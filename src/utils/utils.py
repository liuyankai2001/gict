import pubchempy as pcp
from langchain_core.pydantic_v1 import BaseModel, Field
import json
import os
import dotenv

from langchain_openai import ChatOpenAI

# from config import project_config

from langchain_core.output_parsers import JsonOutputParser
from langchain_core.prompts import ChatPromptTemplate

from src.config import project_config

dotenv.load_dotenv(dotenv_path=project_config.ROOT_DIR / '.env')
os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
os.environ["OPENAI_API_BASE"] = os.getenv("OPENAI_BASE_URL")

class Compound_Info(BaseModel):
    precursor_compound:str = Field(default='前体化合物')
    target_compound:str = Field(default='目标化合物')
    host_cell:str = Field(default='None')

demo_dict = {
    '苯丙氨酸':'phenylalanine',
    '苯甲过氧酸':'benzenecarboperoxic_acid'
}

def parse_user_input(user_input:str)-> dict:
    model = ChatOpenAI(
        model="gpt-4o-mini"
    )
    parser = JsonOutputParser(pydantic_object=Compound_Info)
    prompt = ChatPromptTemplate.from_messages([
        ('system',
         '你是一个合成生物学家，主要任务是根据用户输入，你只需提取用户提到的宿主细胞,前体化合物和目标化合物,并将其转化为英文即可，其他不需要做。另外，如果用户输入的宿主细胞是大肠杆菌，那么就将其翻译为‘ecoli’，如果是酵母，就将其翻译为‘yeast’{format_instructions}'),
        ('human', '{user_input}')
    ])
    chain = prompt | model | parser  # 使用 LangChain 链式调用
    result = chain.invoke({
        "user_input": user_input,
        "format_instructions": parser.get_format_instructions()
    })
    return result

def get_english_name(name):
    return demo_dict[name]

def get_inchikey_by_name(name):
    try:
        # 根据化合物名称搜索（默认返回第一个匹配结果）
        compounds = pcp.get_compounds(name, 'name')
        if compounds:
            return compounds[0].inchikey
        else:
            return None
    except Exception as e:
        print(f"查询出错: {e}")
        return None

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