import requests
import numpy as np
import h5py
import shutil

from six import string_types

base_url = 'http://www.illustris-project.org/api/'

def get(path, params=None, headers=None):
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers=headers)

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically

    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r

class IllustrisHalos(object):
    def __init__(self, api_key, sim_name):
        self.api_key = api_key
        self.sim_name = sim_name
        sims = get(base_url, headers={"api-key": self.api_key})
        sim_names = [s['name'] for s in sims['simulations']]
        sim_index = sim_names.index(self.sim_name)
        self.sim_info = get(sims['simulations'][sim_index]['url'],
                            headers={"api-key": self.api_key})

    def get_snapshot_info(self, snapshot):
        headers = {"api-key": self.api_key}
        if isinstance(snapshot, string_types):
            snap_info = get(self.sim_info["snapshots"]+snapshot, headers=headers)
        else:
            snaps = get(self.sim_info['snapshots'], headers=headers)
            snap_info = get(snaps[snapshot]['url'], headers=headers)
        return snap_info

    def get_halo_data(self, snapshot, subhalos, cutout_request=None):
        """
        Obtain subhalos from an Illustris simulation, adding
        the appropriate header information to the files
        so that yt can read them.

        Parameters
        ----------
        snapshot : integer or string
            Which snapshot to download the subhalos from.
        subhalos : dict, int, or list
            The subhalos to download.
        """
        headers = {"api-key": self.api_key}
        if cutout_request is not None:
            headers.update(cutout_request)

        snap_info = self.get_snapshot_info(snapshot)

        if isinstance(subhalos, dict):
            subs = get(snap_info['subhalos'], params=subhalos, headers=headers)
            num_subs = len(subs['results'])
            subhalos = [subs['results'][i]['id'] for i in range(num_subs)]
        if not isinstance(subhalos, list):
            subhalos = [subhalos]

        subs_info = []
        for subhalo in subhalos:
            subhalo_info = get(snap_info["subhalos"]+"%d" % subhalo, headers=headers)
            subs_info.append(subhalo_info)

        hvals = {}
        hvals["Redshift"] = snap_info["redshift"]
        hvals["Omega0"] = self.sim_info["omega_0"]
        hvals["OmegaLambda"] = self.sim_info["omega_L"]
        hvals["HubbleParam"] = self.sim_info["hubble"]
        hvals["MassTable"] = np.array([0, self.sim_info["mass_dm"], 0, 0, 0, 0])
        hvals["NumFilesPerSnapshot"] = 1
        hvals["Time"] = 1.0/(1.0+snap_info["redshift"])
        hvals["BoxSize"] = self.sim_info["boxsize"]

        filenames = []

        for subhalo_info in subs_info:
            subhalo_url = subhalo_info["cutouts"]["subhalo"]

            halo_filename = get(subhalo_url, headers=headers)

            f = h5py.File(halo_filename, "r+")

            box_width = -1.0

            num_parts = []
            for i in range(6):
                if "PartType%d" % i in f:
                    g = f["/PartType%d" % i]
                    if "ParticleIDs" in g:
                        num_parts.append(g["ParticleIDs"].shape[0])
                        x_min = f["/PartType%d" % i]['Coordinates'][:,0].min()
                        y_min = f["/PartType%d" % i]['Coordinates'][:,1].min()
                        z_min = f["/PartType%d" % i]['Coordinates'][:,2].min()
                        x_width = f["/PartType%d" % i]['Coordinates'][:,0].max() - x_min
                        y_width = f["/PartType%d" % i]['Coordinates'][:,1].max() - y_min
                        z_width = f["/PartType%d" % i]['Coordinates'][:,2].max() - z_min
                        box_width = max(box_width, x_width, y_width, z_width)
                    else:
                        num_parts.append(0)
                else:
                    num_parts.append(0)
            hvals["NumPart_ThisFile"] = np.array(num_parts, dtype="int32")
            hvals["BoxSize"] = 1.1*box_width

            for i in range(6):
                if "PartType%d" % i in f:
                    g = f["/PartType%d" % i]
                    if "Coordinates" in g:
                        f["/PartType%d" % i]['Coordinates'][:,0] -= (subhalo_info['pos_x'] - 0.5*box_width)
                        f["/PartType%d" % i]['Coordinates'][:,1] -= (subhalo_info['pos_y'] - 0.5*box_width)
                        f["/PartType%d" % i]['Coordinates'][:,2] -= (subhalo_info['pos_z'] - 0.5*box_width)

            for k, v in hvals.items():
                f["/Header"].attrs[k] = v

            f.create_group("/Units")

            f["/Units"].attrs["UnitLength_in_cm"] = 3.0856775809623245e+21
            f["/Units"].attrs["UnitMass_in_g"] = 1.98841586e+43
            f["/Units"].attrs["UnitVelocity_in_cm_per_s"] = 1.0e5

            f.flush()
            f.close()

            new_filename = "_".join([self.sim_name, str(snapshot), halo_filename])
            shutil.move(halo_filename, new_filename)

            filenames.append(new_filename)

        return subs_info, filenames
