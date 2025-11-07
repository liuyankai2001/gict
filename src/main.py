import json
import os
import dotenv
import pandas as pd
from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage
from config import project_config
from path_search.get_path_search import get_path_search
from utils.utils import get_english_name, get_inchikey_by_name, get_enzymes_list, parse_user_input, \
    generage_project_file
from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate







def main():
    user_input = "我想要在大肠杆菌中以苯丙氨酸为前体化合物合成苯甲过氧酸，请你给出一段完整的质粒基因序列"
    input_dict = parse_user_input(user_input)


    precursor_compound = input_dict['precursor_compound'].lower().replace(' ', '')
    target_compound = input_dict['target_compound'].lower().replace(' ', '')

    generage_project_file(precursor_compound,target_compound)

    pathways_result_path = project_config.OUTPUT_DIR / str(precursor_compound+'_to_'+target_compound) / 'pathways.csv'

    # print(precursor_compound)
    # print(target_compound)
    # print(pathways_result_path)
    df = None
    if os.path.exists(pathways_result_path):
        df = pd.read_csv(pathways_result_path,nrows=10)
    else:
        get_path_search(str(precursor_compound+'_to_'+target_compound))

    pathway_enzymes = df.loc[0,'pathway_enzymes']
    print(pathway_enzymes)
    enzymes = get_enzymes_list(pathway_enzymes)
    print(enzymes)

if __name__ == '__main__':
    main()