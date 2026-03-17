import os
import shutil
import random
import kim_edn
import subprocess
from add_kimprovenance_entry import add_kimprovenance_entry
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


models_avail = {}
for filename in os.listdir():
    if not os.path.isdir(filename):
        continue
    if filename == "__pycache__":
        continue
    old_extended_id_split = filename.split("_")
    year_and_potential_rank = old_extended_id_split[2]
    potential_rank = year_and_potential_rank[4:]
    species = old_extended_id_split[3]
    if species not in models_avail:
        models_avail[species] = [potential_rank]
    else:
        models_avail[species].append(potential_rank)

recommended_models = {}
with open("recommended_models.yaml") as f:
    recommended_models_by_arity = yaml.load(f, Loader)

for key in recommended_models_by_arity:
    for hyphenated_species in recommended_models_by_arity[key]:
        recommended_models[hyphenated_species.replace("-", "")] = (
            recommended_models_by_arity[key][hyphenated_species]
        )

with open("dois.txt") as f, open("kimid_to_doi.txt", "w") as f_kimid_to_doi:
    for filename in os.listdir():
        if not os.path.isdir(filename):
            continue
        if filename == "__pycache__":
            continue
        
        # Basic stuff -- move, change kimnum, search replace
        new_kimnum = str(random.randint(0, int("9"*12)))
        old_extended_id_split = filename.split("_")
        old_kimnum = old_extended_id_split[-2]
        new_extended_id_split = ["PolyMLP"] + old_extended_id_split[1:-2] + [new_kimnum, "000"]
        new_extended_id = "_".join(new_extended_id_split)
        new_path = f"/home/openkim/models/{new_extended_id}"
        shutil.copytree(filename, new_path)
        subprocess.run(f"sed -i 's/Polymlp/PolyMLP/g' {new_path}/*", shell=True, check=True)
        subprocess.run(f"sed -i 's/{old_kimnum}/{new_kimnum}/g' {new_path}/*", shell=True, check=True)
        new_kimspec_path = f"{new_path}/kimspec.edn"
        kimspec_dict = kim_edn.load(new_kimspec_path)
        doi = next(f)
        kimspec_dict["doi"] = doi.strip()
        kimspec_dict["publication-year"] = "2026"        
        
        # Figure out disclaimer
        species = old_extended_id_split[3]        
        sister_models = sorted(
            sorted(models_avail[species], key=lambda x: x[1]), key=lambda x: len(x)
            )
        year_and_potential_rank = old_extended_id_split[2]
        potential_rank = year_and_potential_rank[4:]
        try:
            recommended_model_and_note = recommended_models[species]
            rec_extended_id_split = recommended_model_and_note["model"].split("_")
            rec_year_and_potential_rank = rec_extended_id_split[2]
            rec_rank = rec_year_and_potential_rank[4:]
            disclaimer = recommended_model_and_note.get("notes")
            if disclaimer is None:
                disclaimer = ""
            else:
                disclaimer += ". "
        except KeyError:
            print(f"No recommended model found for {species}")
            disclaimer = ""
            rec_rank = None
        if len(sister_models) == 1:
            print(f"{species} has only one model available")
        else:
            disclaimer += (
                f"Several models from the PolyMLP repository are available for {species}, "
                "all laying on a pareto front of computational speed vs. RMS accuracy. "
                f"They are labeled {', '.join(sister_models)} in their KIM IDs. "
                "Higher numbers mean slower models"
            )
            if any(["hybrid" in model for model in sister_models]):
                disclaimer += ", and 'hybrid' means that a model is a superposition of polynomial MLPs."
            else:
                disclaimer += "."
            if rec_rank is not None:
                disclaimer += (
                    " The author-recommended model for "
                    f"this system is {rec_rank} "
                    f"({'' if rec_rank == potential_rank else 'not '}this model)."
                )        
        if disclaimer != "":
            kimspec_dict["disclaimer"] = disclaimer 
        kim_edn.dump(kimspec_dict, new_kimspec_path,indent=1)
        add_kimprovenance_entry(new_path, "4ad03136-ed7f-4316-b586-1e94ccceb311", "initial-creation", "")
        print(f"{new_extended_id} {doi.strip()}", file=f_kimid_to_doi)