from Bio import Entrez
from Bio import SeqIO
from Bio.Data import CodonTable
from Bio.Seq import Seq
import time

from src.enzymatic_conversion.codon_usage import ecoli_codon_usage, yeast_coden_usage

# 设置邮箱
Entrez.email = "283293400@qq.com"  # 替换为您的邮箱

# def search_enzyme_CDS(ec_num):
#
#     search_term = f"{ec_num}[EC] AND bacteria[ORGN]"  # 可修改为"Pseudomonas[ORGN]"等特定
#     handle = Entrez.esearch(db="protein", term=search_term, retmax=10)
#     record = Entrez.read(handle)
#     protein_ids = record["IdList"]
#     handle.close()
#     print(f"找到 {len(protein_ids)} 个蛋白ID: {protein_ids}")
#     if not protein_ids:
#         raise ValueError("No protein found for the search term.")
#
#
#     for prot_id in protein_ids:
#         print(f"\n处理蛋白ID: {prot_id}")
#
#         # 步骤2: 链接到核苷酸ID
#         link_handle = Entrez.elink(dbfrom="protein", db="nuccore", id=prot_id, linkname="protein_nuccore")
#         link_record = Entrez.read(link_handle)
#         link_handle.close()
#
#         nuc_ids = []
#         if link_record and link_record[0]["LinkSetDb"]:
#             for link_set in link_record[0]["LinkSetDb"]:
#                 if link_set["DbTo"] == "nuccore":
#                     nuc_ids = [link["Id"] for link in link_set["Link"]]
#
#         print(f"相关核苷酸ID: {nuc_ids}")
#
#         if nuc_ids:
#             nuc_id = nuc_ids[0]  # 取第一个
#
#             # 步骤3: 获取GenBank记录
#             fetch_handle = Entrez.efetch(db="nuccore", id=nuc_id, rettype="gb", retmode="text")
#             gb_record = SeqIO.read(fetch_handle, "genbank")
#             fetch_handle.close()
#
#             # 提取CDS和蛋白序列
#             cds_seq = None
#             protein_seq = None
#             for feature in gb_record.features:
#                 if feature.type == "CDS" and 'translation' in feature.qualifiers:
#                     cds_seq = feature.location.extract(gb_record.seq) # 提取的CDS
#                     protein_seq = Seq(feature.qualifiers['translation'][0]) # 对应的蛋白质序列
#                     break
#
#             if cds_seq is None or protein_seq is None:
#                 continue
#
#             if cds_seq and protein_seq:
#                 print(f"原始CDS: {cds_seq}")
#                 print(f"蛋白序列: {protein_seq}")
#             return cds_seq


def search_enzyme_CDS(ec_num: str, erro_ec,retmax: int = 20):
    """
    根据EC号在细菌中搜索蛋白，返回第一个能成功获取到CDS核苷酸序列的序列。
    如果全部都拿不到，返回 None。
    """
    # search_term = f"{ec_num}[ECNO] AND bacteria[ORGN]"
    # 如果你只想要RefSeq（更可靠），可以改成：
    search_term = f"{ec_num}[ECNO] AND bacteria[ORGN]"
    handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax, sort="relevance")
    record = Entrez.read(handle)
    handle.close()
    protein_ids = record["IdList"]

    if len(protein_ids) == 0:
        search_term = f"{ec_num}[ECNO]"
        handle = Entrez.esearch(db="protein", term=search_term, retmax=retmax, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        protein_ids = record["IdList"]

    print(f"找到 {len(protein_ids)} 个蛋白ID: {protein_ids}",end='')

    if not protein_ids:
        erro_ec.append(ec_num)
        raise ValueError("No protein found for the search term.")

    for prot_id in protein_ids:
        print(f"\n正在尝试蛋白ID: {prot_id}")

        try:
            # 直接取该蛋白精确的CDS核苷酸序列
            fetch_handle = Entrez.efetch(
                db="protein",
                id=prot_id,
                rettype="fasta_cds_na",   # 关键！直接给CDS核苷酸序列
                retmode="text"
            )

            # 使用 parse 而不是 read，防止空返回时出错
            records = list(SeqIO.parse(fetch_handle, "fasta"))
            fetch_handle.close()

            if records:  # 有记录说明成功拿到CDS
                cds_seq = records[0].seq   # 细菌几乎都是单CDS
                print(f"成功！蛋白 {prot_id} 拿到CDS序列，长度 {len(cds_seq)} bp")
                print(f"CDS序列:\n{cds_seq}")
                return cds_seq
            else:
                print(f"蛋白 {prot_id} 没有可用的CDS核苷酸序列，尝试下一个...")

        except Exception as e:
            print(f"蛋白 {prot_id} 请求失败 ({e})，尝试下一个...")

    raise ValueError('所有蛋白都拿不到有效的CDS序列')

def optimize_the_CDS(cds,host_cell):
    """
    根据所选的host_cell，优化原始cds
    :param cds: 酶的原始DNA序列
    :param host_cell: 所选的宿主细胞（如大肠杆菌、酵母细胞）
    :return: 优化后的DNA序列
    """
    original_cds = Seq(cds)
    protein_seq = original_cds.translate(to_stop=True)
    if host_cell=='ecoli':
        codon_usage = ecoli_codon_usage
    elif host_cell=='yeast':
        codon_usage = yeast_coden_usage
    else:
        print('请输入合适的底盘细胞：ecoli or yeast')
        return
    aa_to_best_codon = {}
    for aa, codons in codon_usage.items():
        best_codon = max(codons, key=codons.get)
        aa_to_best_codon[aa] = best_codon

    optimized_cds = 'ATG'  # Start codon
    for aa in protein_seq:
        if aa == '*':  # Skip if stop in protein (but add at end)
            break
        if aa in aa_to_best_codon:
            optimized_cds += aa_to_best_codon[aa]
        else:
            raise ValueError(f"Unknown amino acid in protein sequence: {aa}")
    optimized_cds += aa_to_best_codon['*']  # Add best stop codon

    return optimized_cds

def build_expression_cassette(optimized_cdss, host_cell, promoter_strength, co_expression=True):
    """
    构建表达盒，根据用户选择
    :param optimized_cdss: 优化后的CDS列表 (list of str)
    :param host_cell: 'E. coli' 或 'Yeast'
    :param promoter_strength: 'strong' 或 'weak'
    :param co_expression: True (共表达) 或 False (单独表达)
    :return: 表达盒列表 或 单表达盒字符串
    """
    # 定义元件，根据宿主和强度
    if host_cell == 'ecoli':
        if promoter_strength == 'strong':
            promoter = 'TAATACGACTCACTATAGGG'  # T7 strong
        else:
            promoter = 'TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGC'  # lac weak
        rbs = 'AGGAGGATATCAT'  # E. coli RBS
        terminator = 'AATAAAAAACGCTAGCAGCACTCCGGATGTGCTAGC'  # rrnB T1
    elif host_cell == 'yeast':
        if promoter_strength == 'strong':
            promoter = 'GTTTACAGCTTGTCGGTAAATACCGTTCCAAAGGA'  # PGK1 strong (simplified)
        else:
            promoter = 'GCTTCTTTTTGCTAACGATCAATTTGATCATATG'  # CYC1 weak (simplified)
        rbs = ''  # Yeast不需要RBS，使用 Kozak sequence简化为空或'AACAAA'
        terminator = 'TAAAAAAATTTTTTTTTTT'  # PolyA terminator simplified
    else:
        raise ValueError("Unsupported host")

    if co_expression:
        # 共表达: 共享promoter和terminator，每个CDS前加RBS
        # 启动子 + rbs + cds1 + rbs + cds2 + rbs +...+ rbs + cdsn + 终止子
        cassette = promoter
        for cds in optimized_cdss:
            cassette += rbs + cds
        cassette += terminator
        return cassette
    else:
        # 单独表达: 每个CDS一个完整盒
        cassettes = []
        for cds in optimized_cdss:
            cassette = promoter + rbs + cds + terminator
            cassettes.append(cassette)
        return cassettes

def get_enzymatic_conversion_CDS(enzymes:list[str],host_cell,erro_ec,ec_dict:dict,optimize=False,promoter_strength='strong'):
    origin_cdses = []
    optimize_cds = []
    for enzyme in enzymes:
        print(f"========== 处理酶{enzyme}中... ==========")
        if ec_dict.get(enzyme):
            print(f"{enzyme}已存在于查询字典中")
            ec_cds = ec_dict.get(enzyme)
        else:
            ec_cds = search_enzyme_CDS(enzyme,erro_ec)
            ec_dict[enzyme] = ec_cds
        origin_cdses.append(ec_cds)

    if optimize:
        for ec_cds in origin_cdses:
            opt_cds = optimize_the_CDS(ec_cds,host_cell)
            optimize_cds.append(opt_cds)
    else:
        optimize_cds = origin_cdses
    cassettes = build_expression_cassette(optimize_cds,host_cell,promoter_strength)

    # # 输出结果
    # if co_expression:
    #     print(f"共表达盒: {cassettes}")
    # else:
    #     for i, cass in enumerate(cassettes):
    #         print(f"基因 {i + 1} 表达盒: {cass}")
    #
    return cassettes

if __name__ == '__main__':
    enzymes = ['1.2.3.9', '2.8.3.17']
    cassettes = get_enzymatic_conversion_CDS(enzymes,'ecoli',promoter_strength='strong')
    print(cassettes)