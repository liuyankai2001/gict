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


def is_complete_ec(ec_str):
    parts = ec_str.split('.')
    if len(parts) != 4:
        return False
    # 可选：检查每部分是否为正整数
    try:
        return all(part.isdigit() and int(part) > 0 for part in parts)
    except:
        return False

# def find_first_complete_path(path_list):
#     for path in path_list:
#         ecs = path.split('=>')
#         if all(is_complete_ec(ec) for ec in ecs):
#             return path
#     return None  # 或者 raise ValueError("No complete path found")

def find_ecs_path(path_list):
    for path in path_list:
        ecs = path.split('=>')
        if all(is_complete_ec(ec) for ec in ecs):
            print("酶路径提取成功！酶的路径为：")
            print(path)
            enzymes = get_enzymes_list(path)

    return None  # 或者 raise ValueError("No complete path found")


def get_enzymes_list(enzyme_path):
    return enzyme_path.split('=>')

def generage_project_file(precursor_compound_name,precursor_compound_inchikey,target_compound_name,target_compound_inchikey):
    # 案例
    # precursor_compound_inchikey = 'COLNVLDHVKWLRT'
    # target_compound_inchikey = 'XCRBXWCUXJNEFX'

    input_file = project_config.DEFAULT_PARAMETERS_DIR / 'parameters.txt'

    project_name = str(precursor_compound_name+'_to_'+target_compound_name)
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