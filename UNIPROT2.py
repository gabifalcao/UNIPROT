"""
☆  Origem dos dados:
    UniProt Proteomes (https://www.uniprot.org/proteomes), via API REST pública da UniProt.
    A UniProt calcula e curada os proteomas de referência diretamente a partir de conjuntos de genomas
    anotados e marca explicitamente qual entrada é o “proteoma de referência” para um organismo, com unicidade, e que uma busca de texto no GenBank não consegue fornecer.

☆  Recupera o PROTEOMA DE REFERÊNCIA de um organismo — 
    o conjunto de proteínas curadas, com uma cópia por gene, 
    derivadas de um genoma de referência anotado — ao contrário do genbank_pipeline.py, 
    que recupera um CATÁLOGO de quaisquer registros de proteínas individuais que por acaso tenham sido
    depositados sob o nome desse organismo no GenBank. São coisas realmente diferentes

☆  LIMITAÇÃO:
    Se um organismo não tiver um genoma de referência sequenciado E anotado, 
    não existe proteoma para buscar — de nenhum lugar, por ninguém.

☆  REQUERIMENTO: 
  ✦pip install requests pandas
  ✦quando for botar o organismo que for buscar tem que ser pelo powersheel
        se não, não funciona (de preferencia colocar o nome cientifico) e
        tem que colocar o caminho do arquivo tambem para abrir o codigo,como não é todo mundo
        que sabe disso, exemplo abaixo: 
        cd "C:\Users\GABI\UFMG"(onde o codigo está no meu pc)
        python UNIPROT2.py "Organismo Organismo"  (para rodar e pegar o organismo, escrever de letra maiscula e pode ter espaço)

    """

import os
import logging
from pathlib import Path

import requests

UNIPROT_BASE = "https://rest.uniprot.org"
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

logging.basicConfig(
    filename=OUTPUT_DIR / "pipeline.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",)
logger = logging.getLogger(__name__)

HEADERS = {"User-Agent": f"LBS-proteome-fetcher ({os.getenv('NCBI_EMAIL', 'contact-not-set')})"}


def find_taxon_id(organism_name):
    resp = requests.get(
        f"{UNIPROT_BASE}/taxonomy/search",
        params={"query": organism_name, "format": "json"},
        headers=HEADERS,
        timeout=20,)
    resp.raise_for_status()
    results = resp.json().get("results", [])
    if not results:
        raise ValueError(f"No taxonomy match found for '{organism_name}'")

    
    for r in results:
        if r["scientificName"].lower() == organism_name.lower():#melhor usar o nome cientifico
            return r["taxonId"], r["scientificName"]
    return results[0]["taxonId"], results[0]["scientificName"]


def find_reference_proteome(taxon_id):
    resp = requests.get(
        f"{UNIPROT_BASE}/proteomes/search",
        params={"query": f"organism_id:{taxon_id}", "format": "json"},#procura entrada do proteoma e retorna ele ou nenhum se n existir 
        headers=HEADERS,
        timeout=20,)
    resp.raise_for_status()
    results = resp.json().get("results", [])

    for entry in results:
        if "Reference" in entry.get("proteomeType", ""):
            return entry
    return None

#pega o fasta do proteoma completo diretamente em pedaços por causa do tamanho, para evitar o uso total da memoria
def download_proteome_fasta(upid, out_path):
    url = f"{UNIPROT_BASE}/uniprotkb/stream"
    params = {"query": f"proteome:{upid}", "format": "fasta"}

    with requests.get(url, params=params, headers=HEADERS, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)

    logger.info(f"Saved proteome FASTA to {out_path}")


#colunas para colocar na tabela pq o FASTA não vem com essas informações anotadas
ANNOTATION_FIELDS = [
    "accession", "id", "reviewed", "protein_name", "gene_names",
    "organism_name", "length", "protein_existence", "annotation_score",
    "go_id", "xref_interpro", "xref_pfam", "cc_function",
    "cc_subcellular_location", "ec",
]

#a tabela com a anotação, nome do gene, dominio, função etc, muitas vão estar vazias por falta de estudos e isso é uma informação presente na tabela tambem
def download_annotated_tsv(upid, out_path):
    url = f"{UNIPROT_BASE}/uniprotkb/stream"
    params = {
        "query": f"proteome:{upid}",
        "format": "tsv",
        "fields": ",".join(ANNOTATION_FIELDS),
    }

    with requests.get(url, params=params, headers=HEADERS, stream=True, timeout=120) as resp:
        resp.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)

    logger.info(f"Saved annotated proteome table to {out_path}")


def save_proteome_report(upid, out_path):#proprios metadados do UniProt em um JSON e responde o quanto o proteoma está completo
    import json

    resp = requests.get(
        f"{UNIPROT_BASE}/proteomes/{upid}",
        params={"format": "json"},
        headers=HEADERS,
        timeout=30,)
    resp.raise_for_status()
    data = resp.json()

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)

    busco = (data.get("proteomeCompletenessReport") or {}).get("buscoReport") or {}
    assembly = data.get("genomeAssembly") or {}
    summary = {
        "proteome_id": upid,
        "type": data.get("proteomeType"),
        "strain": data.get("strain"),
        "assembly": assembly.get("assemblyId"),
        "protein_count": data.get("proteinCount"),
        "annotation_score": data.get("annotationScore"),
        "busco_complete_pct": busco.get("score"),
        "busco_fragmented": busco.get("fragmented"),
        "busco_missing": busco.get("missing"),
        "busco_total": busco.get("total"),
        "busco_lineage": busco.get("lineageDb"),
    }
    logger.info(f"Proteome report: {summary}")
    return summary


def fetch_reference_proteome(organism_name):#camimho do arquivo, retorna se não existir tambem

    logger.info(f"Resolving taxon ID for '{organism_name}'")
    taxon_id, matched_name = find_taxon_id(organism_name)
    logger.info(f"Matched '{organism_name}' -> taxon_id={taxon_id} ({matched_name})")

    entry = find_reference_proteome(taxon_id)
    if entry is None:
        logger.warning(f"No reference proteome exists for '{matched_name}' (taxon_id={taxon_id})")
        return None

    upid = entry["id"]
    protein_count = entry.get("proteinCount", "unknown")
    logger.info(f"Found reference proteome {upid} for {matched_name} — {protein_count} proteins")

    safe_name = matched_name.replace(" ", "_")
    base = OUTPUT_DIR / f"{safe_name}_{upid}"
    fasta_path = Path(f"{base}_proteome.fasta")
    tsv_path = Path(f"{base}_annotated.tsv")
    report_path = Path(f"{base}_report.json")

    download_proteome_fasta(upid, fasta_path)
    download_annotated_tsv(upid, tsv_path)
    summary = save_proteome_report(upid, report_path)

    return {
        "fasta": fasta_path,
        "annotated_tsv": tsv_path,
        "report_json": report_path,
        "summary": summary,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Fetch an organism's reference proteome from UniProt."
    )
    parser.add_argument("organism", help='Organism name, e.g. "Escherichia coli"')
    args = parser.parse_args()

    result = fetch_reference_proteome(args.organism)
    if result:
        s = result["summary"]
        print(f"FASTA:           {result['fasta']}")
        print(f"Annotated table: {result['annotated_tsv']}")
        print(f"Proteome report: {result['report_json']}")
        print()
        print(f"Proteome {s['proteome_id']} ({s['type']}), assembly {s['assembly']}")
        print(f"Proteins: {s['protein_count']}   UniProt annotation score: {s['annotation_score']}/5")
        if s["busco_complete_pct"] is not None:
            print(
                f"BUSCO ({s['busco_lineage']}): {s['busco_complete_pct']}% complete, "
                f"{s['busco_fragmented']} fragmented, {s['busco_missing']} missing "
                f"of {s['busco_total']} expected genes"
            )
            if s["busco_complete_pct"] < 90:
                print(
                    "NOTE: this proteome is likely incomplete — report the BUSCO "
                    "score alongside it instead of calling it 'all proteins'."
                )
    else:
        print(
            f"No reference proteome found for '{args.organism}'. "
            f"This organism likely has no annotated reference genome "
            f"assembly — only a GenBank catalog (see genbank_pipeline.py) "
            f"is possible for it."
        )
