import time
from typing import Dict, Tuple

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq
Entrez.email = "283293400@qq.com"  # 必须设置真实邮箱

# 模板配置：{accession: (mcs_start_pos, mcs_end_pos)}
# mcs_start_pos: 插入起始位置（0-based，MCS前）
# mcs_end_pos: 插入结束位置（MCS后，替换整个MCS）
# PLASMID_TEMPLATES: Dict[str, Tuple[int, int]] = {
#     "ecoli": ("KJ412345.1", 860, 880),      # pSEVA231 MCS ~20 bp
#     "yeast": ("MK492243.1", 1235, 1255),    # pYTK MCS ~20 bp
# }


def insert_cassette(cassette_seq: str | Seq,
                   vector_record: SeqRecord,
                   insert_pos: int,
                   cassette_name: str = "My_Expression_Cassette") -> SeqRecord:
    """
    将单个表达盒插入质粒（支持无缝克隆如Gibson/HiFi/In-Fusion）
    返回完整的SeqRecord对象（推荐！可直接保存为GenBank文件，带feature和circular注解）
    """
    # 统一转为Seq对象
    if isinstance(cassette_seq, str):
        cassette_seq = Seq(cassette_seq.upper().replace("\n", "").replace(" ", ""))

    # 插入序列
    new_seq = vector_record.seq[:insert_pos] + cassette_seq + vector_record.seq[insert_pos:]

    # 创建新record（复制所有原始信息）
    new_record = SeqRecord(new_seq)
    new_record.id = f"{vector_record.id}_{cassette_name}"
    new_record.name = vector_record.name + "_new"
    new_record.description = f"{vector_record.description} | Inserted {cassette_name} at {insert_pos+1}"
    new_record.annotations = vector_record.annotations.copy()
    new_record.annotations["topology"] = "circular"  # 保持环状

    # 复制原始features（重要！否则丢失ori、AmpR等信息）
    new_record.features = vector_record.features.copy()

    # 添加表达盒feature（SnapGene里会高亮显示）
    new_feature = SeqFeature(
        FeatureLocation(insert_pos, insert_pos + len(cassette_seq)),
        type="misc_feature",  # 也可以用 "regulatory" 或 "CDS" 根据需要改
        qualifiers={
            "label": cassette_name,
            "note": "User inserted expression cassette (promoter + CDS + terminator)"
        }
    )
    new_record.features.append(new_feature)

    print(f"插入成功！新质粒长度: {len(new_record.seq)} bp")
    return new_record


Entrez.email = "283293400@qq.com"


def fetch_plasmid(accession: str):
    """从NCBI直接下载质粒GenBank文件"""
    print(f"正在下载质粒 {accession}...")
    handle = Entrez.efetch(db="nuccore", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    print(f"下载完成！质粒长度: {len(record.seq)} bp")
    return record


def find_best_insert_position(record) -> int:
    """
    自动寻找最合适的插入位点（返回0-based位置）
    策略优先级：标注MCS → pUC经典序列匹配 → lacZ alpha → 兜底400
    """
    seq = record.seq.upper()

    # 策略1：找标注的MCS（Addgene/SnapGene版本几乎都有）
    for f in record.features:
        qual_str = str(f.qualifiers).lower()
        if any(kw in qual_str for kw in ["mcs", "multiple cloning", "polylinker", "cloning site"]):
            start = int(f.location.start)
            end = int(f.location.end)
            center = start + (end - start) // 2
            print(f"✓ 检测到标注的MCS特征: {start+1}-{end} (1-based)，推荐插入位置: {center + 1}")
            return center

    # 策略2：pUC19/pUC18/pBluescript等经典载体序列匹配（极其精准）
    puc_mcs_patterns = [
        "GAATTCGAGCTCGGTACCCGGG",   # EcoRI到SmaI，长度21bp
        "CGGTACCCGGGATCC",          # KpnI到BamHI
        "CCCGGGATCCTCTAGA",          # SmaI到XbaI（最独特）
    ]
    for pattern in puc_mcs_patterns:
        pos = seq.find(pattern)
        if pos != -1:
            # 经验值：SmaI中心在第一个pattern的 +20 位置（正好是416）
            insert_pos = pos + len(pattern) // 2 + 2  # 微调到SmaI正中间
            print(f"✓ 检测到pUC类经典MCS序列，推荐插入位置（SmaI中心）: {insert_pos + 1} (1-based)")
            return insert_pos

    # 策略3：找lacZ alpha（蓝白筛选载体）
    for f in record.features:
        if f.type == "CDS" and "lacz" in str(f.qualifiers).lower() and len(f.location) < 800:
            start = int(f.location.start)
            suggested = start + 250   # MCS通常在lacZ alpha后半段
            print(f"✓ 检测到lacZ alpha片段，推荐插入位置约: {suggested + 1} (1-based)")
            suggested = min(suggested, len(seq) - 100)  # 防止越界
            return suggested

    # 兜底：绝大多数克隆载体MCS都在400附近
    print("未匹配到任何特征，使用通用兜底位置 400（对pUC19仍安全）")
    return 399  # 0-based


def build_plasmid(cassettes: str, host: str = "ecoli") -> str:
    """
    通用工具：根据宿主构建完整质粒。

    :param cassettes: 表达盒序列（;分隔多个）
    :param host: "ecoli" 或 "yeast"
    :return: 完整DNA序列（大写字符串）
    """
    if host=="ecoli":
        accession = "L09137.2"
    elif host=="yeast":
        accession = "U03451"
    vector = fetch_plasmid(accession)
    insert_pos = find_best_insert_position(vector)
    new_plasmid_record = insert_cassette(cassettes, vector, insert_pos, cassette_name="MyGene")
    return new_plasmid_record

if __name__ == "__main__":
    accession = "L09137.2"  # pUC19（NCBI最原始版）
    # accession = "U03451"         # pRS426（酵母）
    # accession = "AY424263.1"      # pET-28a(+)
    cassette_sequence = "GATCyour_promoter_ATGXXXXXXXXXXXXXXXXXXXXXXXXXTAA_terminator"  # ←基因表达盒序列例子
    vector = fetch_plasmid(accession)
    insert_pos = find_best_insert_position(vector)  # ← 自动找位置
    new_plasmid_record = insert_cassette(cassette_sequence, vector, insert_pos, cassette_name="MyGene")
    # 保存为GenBank（推荐！最专业）
    SeqIO.write(new_plasmid_record, "pUC19_MyGene.gb", "genbank")

    # 如果你只想要纯DNA序列字符串
    full_dna_sequence = str(new_plasmid_record.seq)
    print(full_dna_sequence)  # 或保存fasta也行