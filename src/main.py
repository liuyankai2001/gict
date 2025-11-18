import json
import os
import pickle

import dotenv
import pandas as pd
from Bio import SeqIO
from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage
from config import project_config
from get_plasmid_dna.get_complete_plasmid_dna import build_plasmid
from path_search.get_path_search import get_path_search
from src.enzymatic_conversion.enzymatic_conversion import get_enzymatic_conversion_CDS
# from src.get_plasmid_dna.get_complete_plasmid_dna import get_complete_dna
from utils.utils import get_enzymes_list, generage_project_file, find_ecs_path, is_complete_ec
from utils.langchain_utils import parse_user_input

# from test_class import Test
# test_input = Test()
# host_cell = test_input.host_cell
# precursor_compound_name = test_input.precursor_compound_name
# precursor_compound_inchikey = test_input.precursor_compound_inchikey
# target_compound_name = test_input.target_compound_name
# target_compound_inchikey = test_input.target_compound_inchikey

def main():
    # 用户读入
    user_input = "我想要在大肠杆菌中以苯丙氨酸为前体化合物合成苯甲过氧酸，请你给出一段完整的质粒基因序列"
    if user_input is None:
        user_input = input("请输入您的需求，如：我想要在大肠杆菌中以苯丙氨酸为前体化合物合成苯甲过氧酸，请你给出一段完整的质粒基因序列")

    # 提取用户需求信息
    input_dict = parse_user_input(user_input)

    host_cell = input_dict['host_cell'].lower().replace(' ', '_')
    precursor_compound_name = input_dict['precursor_compound_name'].lower().replace(' ', '_')
    target_compound_name = input_dict['target_compound_name'].lower().replace(' ', '_')
    precursor_compound_inchikey = input_dict['precursor_compound_inchikey']
    target_compound_inchikey = input_dict['target_compound_inchikey']
    # 生成配置文件
    generage_project_file(precursor_compound_name,precursor_compound_inchikey,target_compound_name,target_compound_inchikey)

    # pathways_result_path = project_config.OUTPUT_DIR / str(precursor_compound_name+'_to_'+target_compound_name) / 'pathways.csv'

    # 输出文件路径
    output_path = project_config.OUTPUT_DIR / str(precursor_compound_name+'_to_'+target_compound_name)
    pathways_result_path = output_path / 'pathways.csv'
    output_final_dna_path = output_path / str(host_cell+"_MyGene.gb")
    # pathways_result_path = project_config.OUTPUT_DIR / 'benzenecarboperoxic_acid'/'pathways.csv'

    if os.path.exists(pathways_result_path):
        df = pd.read_csv(pathways_result_path,nrows=30)
    else:
        get_path_search(str(precursor_compound_name+'_to_'+target_compound_name))
        df = pd.read_csv(pathways_result_path, nrows=30)

    pathway_enzymes = df['pathway_enzymes'].tolist()
    # print(pathway_enzymes)
    # enzyme_path = find_ecs_path(pathway_enzymes)
    print("酶路径提取成功！")
    # print(enzyme_path)
    # enzymes = get_enzymes_list(enzyme_path)
    # # 构建基因表达盒
    # print("开始构建基因表达盒...")
    # cassettes = get_enzymatic_conversion_CDS(enzymes, host_cell, promoter_strength='strong')
    erro_ec = []
    ec_dict = {}
    for id,path in enumerate(pathway_enzymes):
        ecs = path.split('=>')
        if all(is_complete_ec(ec) for ec in ecs):
            print(f"尝试第{id+1}条酶路径，酶的路径为：")
            print(path)
            enzymes = get_enzymes_list(path)
            if any(ec in erro_ec for ec in enzymes):
                print("有错误路径存在，跳过")
                continue
            print("开始尝试构建表达盒")
            try:
                # 构建基因表达盒
                cassettes = get_enzymatic_conversion_CDS(enzymes, host_cell, erro_ec,ec_dict,promoter_strength='strong')
            except ValueError as e:
                print(f"出现错误，{e}")

                continue
            else:
                print("成功获取基因表达盒！")
                break
    # print(cassettes)
    # cassettes_path = 'test/data/cassettes.pkl'
    # with open(cassettes_path,'wb') as f:
    #     pickle.dump(cassettes, f)
    # 启动子 + rbs + cds1 + rbs + cds2 + rbs +...+ rbs + cdsn + 终止子
    # p = 'E:\\python_project\\gcit\\src\\test\data\\cassettes.pkl'
    # with open(p, 'rb') as f:
    #     cassettes = pickle.load(f)
    print("基因表达盒构建完毕，开始构建最终DNA序列...")
    final_dna = build_plasmid(cassettes,host_cell)

    SeqIO.write(final_dna, output_final_dna_path, "genbank")
    full_dna_sequence = str(final_dna.seq)
    print(full_dna_sequence)  # 或保存fasta也行

if __name__ == '__main__':
    main()