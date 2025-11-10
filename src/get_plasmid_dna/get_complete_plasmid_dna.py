import time
from typing import Dict, Tuple
from Bio import Entrez, SeqIO
from Bio.Seq import Seq

Entrez.email = "283293400@qq.com"  # 必须设置真实邮箱

# 模板配置：{accession: (mcs_start_pos, mcs_end_pos)}
# mcs_start_pos: 插入起始位置（0-based，MCS前）
# mcs_end_pos: 插入结束位置（MCS后，替换整个MCS）
PLASMID_TEMPLATES: Dict[str, Tuple[int, int]] = {
    "ecoli": ("KJ412345.1", 860, 880),      # pSEVA231 MCS ~20 bp
    "yeast": ("MK492243.1", 1235, 1255),    # pYTK MCS ~20 bp
}


def fetch_plasmid_sequence(accession: str) -> Seq:
    """
    从NCBI拉取质粒序列（FASTA格式，轻量高效）。

    :param accession: NCBI accession
    :return: Bio.Seq.Seq 对象
    """
    try:
        handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        time.sleep(0.34)  # 遵守NCBI速率限制 (~3次/秒)
        return record.seq
    except Exception as e:
        raise ValueError(f"拉取失败 {accession}: {e}")

def insert_cassettes(plasmid_seq: Seq, mcs_start: int, mcs_end: int, cassettes: str) -> Seq:
    """
    将表达盒插入MCS（替换MCS区域）。

    :param plasmid_seq: 模板序列
    :param mcs_start: MCS起始位置
    :param mcs_end: MCS结束位置
    :param cassettes: 表达盒（;分隔）
    :return: 插入后的序列
    """
    # 拆分并拼接多个盒子（添加短间隔如GGGGG，避免框架移位）
    cassette_list = [c.strip().upper() for c in cassettes.split(";") if c.strip()]
    if not cassette_list:
        raise ValueError("无有效表达盒")

    spacer = "GGGGG"  # 5 bp间隔（可选，可移除）
    multi_cassette = spacer.join(cassette_list)

    # 线性化：前 + 插入 + 后
    left = plasmid_seq[:mcs_start]
    right = plasmid_seq[mcs_end:]
    full_seq = left + Seq(multi_cassette) + right

    # 验证：检查长度变化合理
    if len(full_seq) < len(plasmid_seq):
        raise ValueError("插入失败：序列变短")

    return full_seq


def get_complete_dna(cassettes:str):
    """
    将表达盒插入模板质粒，获取完整的质粒序列
    :param cassettes:表达盒
    :return: 完整的质粒DNA序列
    """
    dna = fetch_plasmid_sequence('J01749')
    return dna


def build_plasmid(cassettes: str, host: str = "ecoli") -> str:
    """
    通用工具：根据宿主构建完整质粒。

    :param cassettes: 表达盒序列（;分隔多个）
    :param host: "ecoli" 或 "yeast"
    :return: 完整DNA序列（大写字符串）
    """
    if host.lower() not in PLASMID_TEMPLATES:
        raise ValueError(f"不支持宿主 {host}，可选: ecoli, yeast")

    template_key = host.lower()
    accession, mcs_start, mcs_end = PLASMID_TEMPLATES[template_key]

    # 1. 拉取模板
    print(f"拉取模板: {accession} ({host})...")
    plasmid_seq = fetch_plasmid_sequence(accession)

    # 2. 插入表达盒
    print(f"插入 {len(cassettes.split(';'))} 个表达盒...")
    final_seq = insert_cassettes(plasmid_seq, mcs_start, mcs_end, cassettes)

    # 3. 返回字符串（环状，假设序列已闭合）
    print(f"构建完成！长度: {len(final_seq)} bp")
    return str(final_seq).upper()


if __name__ == "__main__":
    # 示例表达盒（替换为你的真实序列）
    example_cassettes = (
        "TTATAATGTATAATGCTAGCTTCTGCGCTCGCTCGCTACTAGAGGCTAGCGGAATTCGAGCTCGGTACCCGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTTGGCGAATTGCGGCCGCGTCTAGAAGGCCGGCGCCATTGATGTAGCGGTAGCG"  # 盒1: 启动子+CDS+终止子
        ";"
        "T7_promoter_ATG_GFP_TAA_polyA"  # 盒2: 简化
    )

    # 大肠杆菌构建
    ecoli_plasmid = build_plasmid(example_cassettes, host="ecoli")
    print("E. coli 质粒片段:", ecoli_plasmid[:200], "...")

    # 酵母构建
    yeast_plasmid = build_plasmid("你的酵母表达盒", host="yeast")
    print("Yeast 质粒片段:", yeast_plasmid[:200], "...")

    # 保存到文件（可选）
    # with open("ecoli_plasmid.fasta", "w") as f:
    #     f.write(f">pSEVA231_inserted\n{ecoli_plasmid}")