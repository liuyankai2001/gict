from dataclasses import dataclass

@dataclass
class Test:
    host_cell:str = "ecoli"
    precursor_compound_name:str = "phenylalanine"
    target_compound_name:str = "benzoyl_peroxide"
    precursor_compound_inchikey:str = "COLNVLDHVKWLRT"
    target_compound_inchikey:str = "XCRBXWCUXJNEFX"
