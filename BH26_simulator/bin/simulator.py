import os
import sys
import array

sys.path.append("../../code")
from cellDynamics import CellDynamicsMosquitoBH26Delay

class Simulator:

    def genotype_string(self, genotype_num):
        if genotype_num == 0:
            return "ww"
        elif genotype_num == 1:
            return "wc"
        elif genotype_num == 2:
            return "wr"
        elif genotype_num == 3:
            return "cc"
        elif genotype_num == 4:
            return "cr"
        elif genotype_num == 5:
            return "rr"
        return 'None'

    def reduction_matrix(self, params):
        matrix = []
        for gM in range(6):
            matrix.append([])
            for gF in range(6):
                matrix[gM].append(params[gM * 2] * params[1 + gF * 2])
        return matrix
                
    def reduction_params(self, matrix):
        params = []
        for genotype in range(1, 6):
            params.append(matrix[genotype][0]) # male genotype with female ww
            params.append(matrix[0][genotype]) # male ww with female genotype
        return params

    def __init__(self, fn):
        with open(fn, 'r') as f:
            data = f.readlines()
        self.values = {}
        for key in ["num_species", "delay", "death_rate", "competition", "emergence_rate", "activity", "reduction_parameters", "hybridisation", "sex_ratio", "female_bias", "m_w", "m_c", "small_value", "initial_populations", "qm", "day_of_introduction", "species_number_of_introduction", "genotype_of_introduction", "sex_of_introduction", "quantity_of_introduction", "simulation_days"]:
            self.values[key] = None
        
        for line in data:
            if not line.strip(): continue
            line = line.strip()
            if line.startswith("#"): continue
            key, line = line.split(",", 1)
            key = key.strip()

            num_species = self.values["num_species"]
            if num_species == None: num_species = 3
            if key == "num_species" or key == "delay" or key == "day_of_introduction" or key == "simulation_days" or key == "quantity_of_introduction":
                self.values[key] = int(line.split(",")[0])
            elif key == "species_number_of_introduction":
                self.values[key] = int(line.split(",")[0]) - 1
            elif key == "sex_of_introduction":
                line = line.split(",")[0].strip().lower()
                if line == "male":
                    self.values[key] = 0
                elif line == "female":
                    self.values[key] = 1
                else:
                    sys.stderr.write("Error: sex_of_introduction must be either 'Male' or 'Female'\n")
                    sys.exit(1)
            elif key == "genotype_of_introduction":
                line = line.split(",")[0].strip().lower()
                if line == "ww":
                    self.values[key] = 0
                elif line == "wc":
                    self.values[key] = 1
                elif line == "wr":
                    self.values[key] = 2
                elif line == "cc":
                    self.values[key] = 3
                elif line == "cr":
                    self.values[key] = 4
                elif line == "rr":
                    self.values[key] = 5
                else:
                    sys.stderr.write("Error: genotype_of_introduction must be ww, wc, wr, cc, cr or rr\n")
                    sys.exit(1)
            elif key == "death_rate":
                s = list(map(float, line.split(",")[:2 * 6 * num_species]))
                male = s[0::2]
                female = s[1::2]
                self.values["death_rate"] = [[male[i::6] for i in range(6)], [female[i::6] for i in range(6)]]
            elif key == "competition":
                s = list(map(float, line.split(",")[:num_species * num_species]))
                self.values["competition"] = [s[i::3] for i in range(3)]
            elif key == "emergence_rate":
                self.values["emergence_rate"] = list(map(float, line.split(",")[: num_species]))
            elif key == "activity":
                s = list(map(float, line.split(",")[:num_species * num_species]))
                self.values["activity"] = [s[i::3] for i in range(3)]
            elif key == "reduction_parameters":
                s = list(map(float, line.split(",")[:2 * (6 - 1)]))
                self.values[key] = [1, 1] + s
            elif key == "hybridisation":
                s = list(map(float, line.split(",")[:num_species * num_species * num_species]))
                s = [s[i::3] for i in range(3)]
                self.values["hybridisation"] = [[x[i::3] for i in range(3)] for x in s]
            elif key == "sex_ratio" or key == "female_bias" or key == "m_w" or key == "m_c" or key == "small_value":
                self.values[key] = float(line.split(",")[0])
            elif key == "initial_populations":
                self.values[key] = list(map(float, line.split(",")[:2 * num_species]))
            elif key == "qm":
                self.values[key] = list(map(float, line.split(",")[:num_species]))
        self.cell = CellDynamicsMosquitoBH26Delay(num_species = self.values["num_species"],
                                                  delay = self.values["delay"],
                                                  current_index = 0,
                                                  death_rate = self.values["death_rate"],
                                                  competition = self.values["competition"],
                                                  emergence_rate = self.values["emergence_rate"],
                                                  activity = self.values["activity"],
                                                  reduction = self.reduction_matrix(self.values["reduction_parameters"]),
                                                  hybridisation = self.values["hybridisation"],
                                                  sex_ratio = self.values["sex_ratio"],
                                                  female_bias = self.values["female_bias"],
                                                  m_w = self.values["m_w"],
                                                  m_c = self.values["m_c"])
        self.cell.setSmallValue(self.values["small_value"])

        self.populations = []
        

    def getValues(self):
        return self.values

    def getPopulations(self):
        return self.populations

    def getCurrentAdults(self, pops_and_params):
        num_species = self.values["num_species"]
        num_genotypes = 6
        num_sexes = 2
        current_adults = [0 for i in range(num_species * num_genotypes * num_sexes)]
        for sex in range(num_sexes):
            for genotype in range(num_genotypes):
                for species in range(num_species):
                    for d in [self.cell.getCurrentIndex()]:
                        indi = species + genotype * num_species + sex * num_species * num_genotypes + d * num_species * num_genotypes * num_sexes
                        indj = species + genotype * num_species + sex * num_species * num_genotypes
                        current_adults[indj] = pops_and_params[indi]
        return current_adults


    def go(self):
        num_pops = self.cell.getNumberOfPopulations()
        num_species = self.values["num_species"]
        num_genotypes = 6
        num_sexes = 2
        delay = self.cell.getDelay()

        # initialise adult populations
        initial_condition = [0 for i in range(num_pops)]
        for sex in range(num_sexes):
            for genotype in range(1): # only ww initially
                for species in range(num_species):
                    for d in range(delay + 1):
                        indi = species + genotype * num_species + sex * num_species * num_genotypes + d * num_species * num_genotypes * num_sexes
                        indj = sex + species * num_sexes
                        initial_condition[indi] = self.values["initial_populations"][indj]
        initial_condition += self.values["qm"]
        pap = array.array('f', initial_condition)

        self.populations = []

        # burn-in
        for timestep in range(self.values["day_of_introduction"]):
            self.populations.append(self.getCurrentAdults(pap))
            self.cell.evolve(1.0, pap)
            self.cell.incrementCurrentIndex()

        # introduce
        indi = self.values["species_number_of_introduction"] + self.values["genotype_of_introduction"] * num_species + self.values["sex_of_introduction"] * num_species * num_genotypes + self.cell.getCurrentIndex() * num_species * num_genotypes * num_sexes
        pap[indi] += self.values["quantity_of_introduction"]

        # evolve
        for timestep in range(self.values["day_of_introduction"], self.values["simulation_days"]):
            self.populations.append(self.getCurrentAdults(pap))
            self.cell.evolve(1.0, pap)
            self.cell.incrementCurrentIndex()

        
        
    def output(self, fn):
        num_species = self.values["num_species"]
        with open(fn, 'w') as f:
            f.write("# Output produced by BH26Delay simulator with the following parameters\n")
            f.write("num_species, " + str(num_species) + "\n")
            f.write("delay, " + str(self.cell.getDelay()) + "\n")

            comment = "# Death rates"
            vals = "death_rate"
            dr = self.cell.getDeathRate()
            for sp in range(num_species):
                for g in range(6):
                    for sex in range(2):
                        comment += ", " + ("M" if sex == 0 else "F") + self.genotype_string(g) + "_" + str(sp + 1)
                        vals += ", " + str(dr[sex][g][sp])
            f.write(comment + "\n")
            f.write(vals + "\n")

            comment = "# Competition"
            vals = "competition"
            c = self.cell.getCompetition()
            for s1 in range(num_species):
                for s2 in range(num_species):
                    comment += ", " + str(s1 + 1) + "v" + str(s2 + 1)
                    vals += ", " + str(c[s2][s1])
            f.write(comment + "\n")
            f.write(vals + "\n")
            
            comment = "# Emergence rate"
            vals = "emergence_rate"
            c = self.cell.getEmergenceRate()
            for s1 in range(num_species):
                comment += ", " + str(s1 + 1)
                vals += ", " + str(c[s1])
            f.write(comment + "\n")
            f.write(vals + "\n")

            comment = "# activity between species F=Female M=Male"
            vals = "activity"
            c = self.cell.getActivity()
            for sM in range(num_species):
                for sF in range(num_species):
                    comment += ", F" + str(sF + 1) + "_M" + str(sM + 1)
                    vals += ", " + str(c[sF][sM])
            f.write(comment + "\n")
            f.write(vals + "\n")

            comment = "# reduction parameters (F=Female M=Male)"
            for g in range(1, 6):
                for sex in range(2):
                    comment += ", " + ("M" if sex == 0 else "F") + self.genotype_string(g)
            f.write(comment + "\n")
            f.write("reduction_parameters, " + ",".join(map(str, self.reduction_params(self.cell.getReduction()))) + "\n")

            comment = "# hybridisation F=Female M=Male o=offspring"
            vals = "hybridisation"
            c = self.cell.getHybridisation()
            for so in range(num_species):
                for sF in range(num_species):
                    for sM in range(num_species):
                        comment += ", M" + str(sM) + "_F" + str(sF) + "_o" + str(so)
                        vals += ", " + str(c[sM][sF][so])
            f.write(comment + "\n")
            f.write(vals + "\n")

            f.write("sex_ratio, " + str(self.cell.getSexRatio()) + "\n")
            f.write("female_bias, " + str(self.cell.getFemaleBias()) + "\n")
            f.write("m_w, " + str(self.values["m_w"]) + "\n")
            f.write("m_c, " + str(self.values["m_c"]) + "\n")
            f.write("small_value, " + str(self.cell.getSmallValue()) + "\n")

            f.write("# initial wild-types F=Female M=Male, M1, F1, M2, F2, M3, F3 ...\n")
            f.write("initial_populations, " + ",".join(map(str, self.values["initial_populations"])) + "\n")
            f.write("# qm values (related to carrying capacity), QM1, QM2, QM3 ...\n")
            f.write("qm, " + ",".join(map(str, self.values["qm"])) + "\n")
            f.write("day_of_introduction, " + str(self.values["day_of_introduction"]) + "\n")
            f.write("genotype_of_introduction, " + self.genotype_string(self.values["genotype_of_introduction"]) + "\n")
            f.write("species_number_of_introduction, " + str(self.values["species_number_of_introduction"] + 1) + "\n")
            f.write("sex_of_introduction, " + ("Male" if (self.values["sex_of_introduction"] == 0) else "Female") + "\n")
            f.write("quantity_of_introduction, " + str(self.values["quantity_of_introduction"]) + "\n")
            f.write("simulation_days, " + str(self.values["simulation_days"]) + "\n")

            f.write("#\n# Results.  First number is species.  Next letters indicate genotype.  Final is Male or Female\n")
            header = "day "
            for species in range(num_species):
                for genotype in range(6):
                    for sex in ['M', 'F']:
                        header += ", " + str(species + 1) + self.genotype_string(genotype) + sex
            f.write(header + "\n")
            t = 0
            for p in self.populations:
                result = str(t) 
                for species in range(num_species):
                    for genotype in range(6):
                        for sex in range(2):
                            indj = species + genotype * num_species + sex * num_species * 6
                            result += ", " + str(p[indj])
                f.write(result + "\n")
                t += 1
                            
            
