import pubchempy as pcp
from langchain_core.pydantic_v1 import BaseModel, Field
import json
import os
import dotenv
import pandas as pd
from langchain_openai import ChatOpenAI

# from config import project_config

from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate

from src.config import project_config




def get_enzymes_list(str:str):
    enzymes = str.split('=>')
    return enzymes

def generage_project_file(precursor_compound,target_compound):
    # 案例
    precursor_compound_inchikey = 'COLNVLDHVKWLRT'
    target_compound_inchikey = 'XCRBXWCUXJNEFX'

    input_file = project_config.DEFAULT_PARAMETERS_DIR / 'parameters.txt'

    project_name = str(precursor_compound+'_to_'+target_compound)
    if not os.path.exists(project_config.PROJECT_DIR/ project_name):
        os.makedirs(project_config.PROJECT_DIR/ project_name)

    output_file = project_config.PROJECT_DIR/ project_name / 'parameters.txt'
    print(output_file)

    with open(input_file, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    updated_lines = []
    for line in lines:
        stripped = line.strip()
        # 跳过注释行和空行进行匹配（但保留原格式）
        if stripped.startswith('source|'):
            updated_lines.append(f'source|{precursor_compound_inchikey}\n')
        elif stripped.startswith('target|'):
            updated_lines.append(f'target|{target_compound_inchikey}\n')
        else:
            updated_lines.append(line)  # 保持原有格式（包括换行）

    with open(output_file, 'w', encoding='utf-8') as f:
        f.writelines(updated_lines)


if __name__ == '__main__':
    pass