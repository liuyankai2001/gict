from Bio import Entrez, SeqIO
from Bio.Seq import Seq


def fetch_plasmid_sequence(accession):
    """
    从NCBI获取质粒模板序列。
    :param accession: NCBI accession number (e.g., 'J01749' for pBR322)
    :return: DNA序列字符串
    """
    Entrez.email = "283293400@qq.com"  # 替换为你的邮箱
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        return str(record.seq)
    except Exception as e:
        print(f"Error fetching sequence: {e}")
        return None  # 或手动粘贴序列


def insert_cassette_into_plasmid(plasmid_seq, cassette_seq, insert_site_motif='GAATTC'):
    """
    将表达盒插入质粒模板的特定位点。
    :param plasmid_seq: 质粒模板序列 (str)
    :param cassette_seq: 表达盒序列 (str)
    :param insert_site_motif: 插入位点序列 (e.g., 'GAATTC' for EcoRI)
    :return: 插入后的完整质粒序列 (str)
    """
    position = plasmid_seq.find(insert_site_motif)
    if position == -1:
        raise ValueError(f"Insertion site '{insert_site_motif}' not found in plasmid sequence.")

    # 模拟切开位点后插入 (假设切后插入, 无疤痕; 实际中调整端)
    complete_seq = plasmid_seq[:position + len(insert_site_motif) // 2] + cassette_seq + plasmid_seq[position + len(
        insert_site_motif) // 2:]

    return complete_seq

def get_complete_dna(cassettes:list[str]):
    pass


if __name__ == '__main__':
    print(fetch_plasmid_sequence('J01749'))