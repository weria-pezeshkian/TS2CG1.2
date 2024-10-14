from pathlib import Path
import argparse
import os
import shutil
import numpy as np

class PointUpdaterClass:
    def __init__(self, path="point/"):
        self.path = Path(path)
        self.box = None
        self.input_file = None
        self.lipid_domains = None
        self.inclusions = None
        self.exclusions = None
        self.k = 1

        try:
            self.inner = self._read_bm(self.path / "InnerBM.dat")
            self.outer = self._read_bm(self.path / "OuterBM.dat")

            inclusions_path = self.path / "IncData.dat"
            if inclusions_path.is_file():
                self.inclusions = self._read_data(inclusions_path, coords=True)

            exclusions_path = self.path / "ExcData.dat"
            if exclusions_path.is_file():
                self.exclusions = self._read_data(exclusions_path, coords=False)
        except Exception as e:
            raise ValueError(f"Error initializing PointUpdater: {e}")

    def _read_bm(self, file):
        # Check if there is box definition in the file
        skiprows = 3
        if not self.box:
            with file.open("r") as f:
                for line in f:
                    if "Box" in line:
                        self.box = line.strip()
                        skiprows = 4
                        break

        # Read file as dictionary
        data = np.loadtxt(file, skiprows=skiprows)
        names = "id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2".split()
        return {name: data[:, i] for i, name in enumerate(names)}

    def _read_data(self, file, coords=True):
        data = np.loadtxt(file, skiprows=2)
        names = "id typeid pointid lx ly lz".split() if coords else "id pointid r".split()
        return {name: data[:, i] for i, name in enumerate(names)}

    def update_domains(self, layer = "both", method = None, assign = None):
        domains = {}
        if layer == "both":
            domains = {"inner": self.inner["domain_id"], "outer": self.outer["domain_id"]}
        elif layer in ["inner", "outer"]:
            domains = {layer: getattr(self, layer)["domain_id"]}
        else:
            raise ValueError(f"Invalid layer: {layer}. Choose 'both', 'inner', or 'outer'.")

        if method:
            for key in domains:
                new_domain = method()
                if len(new_domain) != len(domains[key]):
                    raise ValueError(f"Method returned incorrect length: {len(new_domain)}, expected {len(domains[key])}")
                domains[key] = new_domain
        elif assign is not None:
            for key in domains:
                if len(assign) != len(domains[key]):
                    raise ValueError(f"Assigned array has incorrect length: {len(assign)}, expected {len(domains[key])}")
                domains[key] = np.asarray(assign)

        for key, value in domains.items():
            getattr(self, key.lower())["domain_id"] = value

    def assign_by_curvature(self, input_file, location = "both", use_c12 = False):
        self.input_file = Path(input_file)

        self.lipid_domains = self._parse_input_file(use_c12)

        locations = [self.inner, self.outer] if location == "both" else [getattr(self, location)]

        for loc in locations:
            N = len(loc["id"])
            randomizer = np.random.permutation(N)

            for index in randomizer:
                domain = self._assign_domain(loc, self.lipid_domains, index, use_c12)
                loc["domain_id"][index] = domain

    def _parse_input_file(self, use_c12):
        lipid_domains = {}

        if use_c12:
            columns = "lipid percentage APL c1 c2".split()
        else:
            columns = "lipid percentage APL c0".split()

        ndx = 0
        with self.input_file.open("r") as f:
            for line in f:
                if not line.startswith(";"):
                    parts = line.strip().split()
                    lipid_domains[ndx] = {column: value for column, value in zip(columns, parts)}
                    ndx += 1
        return lipid_domains

    def _assign_domain(self, loc, lipid_domains, index, use_c12):
        Cs = np.array([loc["C1"][index], loc["C2"][index]])

        if use_c12:
            for key, value in lipid_domains.items():
                c12 = np.array([float(value.get(c)) for c in ['c1', 'c2']])
                distances = {key: np.linalg.norm(Cs - c12) for key, value in lipid_domains.items()}
            domain = min(distances, key=distances.get)
        else:
            H = 1/2 * Cs.sum()
            for key, value in lipid_domains.items():
                c0 = float(value['c0'])
                probabilities = {key: np.exp(-self.k * (2*H - c0)**2)}
            total_prob = sum(probabilities.values())
            probabilities = {key: value / total_prob for key, value in probabilities.items()}
            domain = np.random.choice(list(probabilities.keys()), p=list(probabilities.values()))

        return int(domain)

    def write_folder(self):
        backup = self._create_backup(self.path)
        self.path.mkdir(parents=True, exist_ok=True)
        print(f"Backup created as {backup}")

        for name, data in [("InnerBM.dat", self.inner), ("OuterBM.dat", self.outer)]:
            self._write_bm_file(name, data)

        if self.inclusions:
            self._write_inc_file("IncData.dat", self.inclusions)

        if self.exclusions:
            self._write_inc_file("ExcData.dat", self.exclusions)

    def _create_backup(self, path):
        backup = path
        while backup.exists():
            backup = Path(f"#{backup}#")
        os.rename(path, backup)
        return backup

    def _write_bm_file(self, name, data):
        file_path = self.path / name

        all_data = np.column_stack([data[key] for key in data])
        header = f"< Point NoPoints {len(data['id'])} >\n< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2 >\n< {name[:5]} >"
        if "Outer" in name:
            header = f"{self.box}\n{header}"
        np.savetxt(file_path, np.round(all_data, 3), header=header, comments='', fmt="".join(['%10d', '%5d'] + ['%10.3f']*4 + ['%8.3f']*11))

    def _write_inc_file(self, name, data):
        file_path = self.path / name

        all_data = np.column_stack([data[key] for key in data])
        header = f"< Inclusion NoInc {len(data['id'])} >\n< id typeid pointid lx ly lz >"
        np.savetxt(file_path, np.round(all_data, 3), header=header, comments='', fmt="".join(['%12d']*3 + ['%8.3f']*3))

    def write_input_str(self, output_str, input_str):

        with output_str.open("w") as out_f, input_str.open("r") as in_f:
            in_lipids_list = False
            for line in in_f:

                if not in_lipids_list:
                    out_f.write(line)

                if "[Lipids List]" in line:
                    for domain, lipid in self.lipid_domains.items():
                        out_f.write(f"Domain {domain}\n{lipid['lipid']} 1 1 {lipid['APL']}\nEnd\n")
                    in_lipids_list = True
                    continue
                elif "End" in line and in_lipids_list:
                    in_lipids_list = False

        # Backup if we are overwriting the str file
        if output_str == input_str:
            backup = self._create_backup(input_str)
            print(f"Backup created as {backup}")

    @classmethod
    def from_args(cls, args):
        parser = argparse.ArgumentParser(description="Update point folder based on curvature preferences")
        parser.add_argument('-i', '--input', type=str, default="domain.txt", help="Path to the domain input file")
        parser.add_argument('-c', '--c12', action='store_true', help="Use c12 approach instead of c0")
        parser.add_argument('-p', '--path', default="point/", help="Path to the point folder")
        parser.add_argument('-l', '--location', default='both', choices=['both', 'upper', 'lower'], help="Choose which monolayer to alter")
        parser.add_argument('-ni', '--new_TS2CG', default=None, help="Path to write a new TS2CG input file", type=Path)
        parser.add_argument('-oi', '--old_TS2CG', default=None, help="Path to the old TS2CG input.str", type=Path)
        parser.add_argument('-I', '--ignore_lipid_number', action='store_true', help="Ignore the number of lipids in the input file")

        parsed_args = parser.parse_args(args)

        updater = cls(parsed_args.path)
        updater.assign_by_curvature(parsed_args.input, parsed_args.location, parsed_args.c12)

        if parsed_args.new_TS2CG:
            updater.write_input_str(parsed_args.new_TS2CG, parsed_args.old_TS2CG)

        updater.write_folder()
        return updater
