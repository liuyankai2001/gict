import requests
from langchain_core.pydantic_v1 import BaseModel, Field
import os
import dotenv
import json
from langchain_core.prompts import ChatPromptTemplate,MessagesPlaceholder
from langchain_openai import ChatOpenAI
from langchain_core.tools import tool
from langchain.agents import create_tool_calling_agent,AgentExecutor
from config import project_config
from langchain_core.output_parsers import JsonOutputParser
from pydantic import BaseModel

class Compound_Info(BaseModel):
    host_cell:str = Field(description='底盘细胞的英文，ecoli or yeast')
    precursor_compound_name:str = Field(description='前体化合物的英文名')
    precursor_compound_inchikey: str = Field(description='前体化合物的inchikey')
    target_compound_name: str = Field(description='目标化合物的英文名')
    target_compound_inchikey: str = Field(description='目标化合物的inchikey')

dotenv.load_dotenv(dotenv_path=project_config.ROOT_DIR / '.env')
os.environ["OPENAI_API_KEY"] = os.getenv("OPENAI_API_KEY")
os.environ["OPENAI_API_BASE"] = os.getenv("OPENAI_BASE_URL")


system_prompt = (
    "你是一个合成生物学家，负责从用户输入中精确提取并转换以下5个字段并忽略用户其他请求：\n\n"
    "1. host_cell: 如果提到大肠杆菌，必须输出 'ecoli'；如果提到酵母，必须输出 'yeast'。\n"
    "2. precursor_compound_inchikey: 必须是前体化合物的 InChIKey 的前14个字符（例如：COLNVLDHVKWLRT）。\n"
    "3. target_compound_inchikey: 必须是目标化合物的 InChIKey 的前14个字符（例如：OMPJBNCRMGITSC）。\n"
    "4. precursor_compound_name: 必须是前体化合物的英文名。\n"
    "5. target_compound_name: 必须是目体化合物的英文名。\n\n"
    "重要规则：\n"
    "- 你必须先将提取出的中文化合物名翻译为英文,然后后传入 `get_inchikey_by_name_api` 工具，。\n"
    "- 工具返回的是完整 InChIKey（如 'COLNVLDHVKWLRT-QMMMGPOBSA-N'），你只需取**前14个字符**。\n"
    "- 在最终 JSON 中，**precursor_compound_inchikey 和 target_compound_inchikey 字段绝对不能是英文名称**，必须是 InChIKey 前缀。\n"
    "- 在最终 JSON 中，precursor_compound_name和target_compound_name必须是对应的前体化合物和目体化合物的英文名。\n"
    "- 最终输出必须严格符合以下 JSON 格式，不要任何额外文字、解释或 Markdown：\n"
    "{format_instructions}"
)

# @tool(description='该函数是一个查询函数，根据所输入的化合物英文名，查询ncbi数据库，得到InChIKey字符')
@tool
def get_inchikey_by_name_api(name:str):
    """
    该函数是一个查询函数，根据所输入的化合物英文名，查询ncbi数据库，得到InChIKey字符
    :param name: sr类型，化合物英文名称
    :return: InChIKey字符串
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/InChIKey/JSON"
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if "PropertyTable" in data:
            return data["PropertyTable"]["Properties"][0]["InChIKey"]
    return None

chat_template = ChatPromptTemplate.from_messages([
    # ('system','你是一个合成生物学家，现在做一个命名实体提取任务。'
    #           '根据用户所输入的一段内容，提取底盘细胞、前体化合物和目标化合物即可，忽略其他要求。'
    #           '将提取底盘细胞、前体化合物和目标化合物翻译成英文，并使用查询工具查询前体化合物和目标化合物的inchikey，仅输出 InChIKey 的前14个字符，不要任何其他文字、标点、空格、说明或格式。'
    #           '这里，如果底盘细胞是大肠杆菌，就翻译为ecoli，如果是酵母，就翻译为yeast。'
    #           '{format_instructions}'),
    ('system',system_prompt),
    ('human','{input}'),
    MessagesPlaceholder('agent_scratchpad')
])

def parse_user_input(user_input):
    llm = ChatOpenAI(model="gpt-4o-mini")
    tools = [get_inchikey_by_name_api]
    parser = JsonOutputParser(pydantic_object=Compound_Info)
    agent = create_tool_calling_agent(llm, prompt=chat_template, tools=tools)
    agent_executor = AgentExecutor(agent=agent, tools=tools, verbose=True)
    res = agent_executor.invoke(
        {"input": user_input,
         "format_instructions": parser.get_format_instructions()})
    output_str = res["output"]
    ans = parser.invoke(output_str)
    return ans

if __name__ == '__main__':
    user_input = '我想要在大肠杆菌中以苯丙氨酸为前体化合物合成过氧化苯甲酰，请你给出一段完整的质粒基因序列'
    input_dict = parse_user_input(user_input)
    print(input_dict)
    print(type(input_dict))