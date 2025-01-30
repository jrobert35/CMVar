import datetime
import os
import pandas as pd
import math
import argparse
from jinja2 import Environment, FileSystemLoader

# Prise en compte des arguments
parser = argparse.ArgumentParser(description="Script de création du rapport HTML de l'analyse CMV")
parser.add_argument("--template_dir", "-t", required=True, help="Chemin vers le repertoire contenant le template HTML")
parser.add_argument("--bed", "-b", required=True, help="Chemin vers le fichier BED")
parser.add_argument("--vcf", "-v", required=True, help="Chemin vers le fichier VCF (sortie de snpeff)")
parser.add_argument("--flagstat", "-f", required=True, help="Chemin vers le fichier flagstat")
parser.add_argument("--coverage", "-c", required=True, help="Chemin vers le fichier de couverture")
parser.add_argument("--output", '-o', required=True, help="Chemin vers le fichier de sortie")
args = parser.parse_args()

# Déclaration des fichiers
vcf_file = args.vcf
flagstat = args.flagstat
bed = args.bed
coverage = args.coverage
output_file = args.output

# Définition de fonction
def get_sentence(df_cov, ul):
    phrases = []
    for index, row in df_cov.iterrows():
        if row[1] == ul:
            start = row['start']
            end = row['end']
            phrase = phrase = "Région non couverte de {} à {}".format(start, end)
            phrases.append(phrase)
    return "<br>".join(phrases)

def check_mutation_presence(row, ul):
    mutations = row["Acides Aminés"].split("/")
    for mutation in mutations:
        if mutation in mutations_ul.get(ul, []):
            return "Oui"
    return ""

# Traitement des fichiers
# Récupération du nom de l'échantillon
ech = os.path.splitext(os.path.basename(vcf_file))[0].replace("_snpeff", "")

# Récupération du nombre de reads + Récupération du pourcentage de reads mappés
with open(flagstat, 'r') as fichier:
    premiere_ligne = fichier.readline()
    mots = premiere_ligne.split()
    if mots:
        sequences = int(mots[0])
	
with open(flagstat, 'r') as fichier:
	lignes = fichier.readlines()
	for ligne in lignes:
		if "mapped" in ligne:
			mots = ligne.split()
			mapped_percentage = mots[4].replace('(', '')
			break

# Récupération des zones non couvertes
if os.path.getsize(coverage) > 0:
    df_bed = pd.read_csv(bed, delimiter='\t', header=None)
    df_cov = pd.read_csv(coverage, delimiter='\t', header=None)

    for index, row in df_bed.iterrows():
        if row[5] == "-" and row[3] in df_cov[1].values:
            df_cov.loc[df_cov[1] == row[3], "start"] = (row[2] - df_cov.loc[df_cov[1] == row[3], 4]) / 3
            df_cov.loc[df_cov[1] == row[3], "end"] = (row[2] - df_cov.loc[df_cov[1] == row[3], 3]) / 3
        elif row[5] == "+" and row[3] in df_cov[1].values:
            df_cov.loc[df_cov[1] == row[3], "start"] = (df_cov.loc[df_cov[1] == row[3], 3] - row[1]) / 3
            df_cov.loc[df_cov[1] == row[3], "end"] = (df_cov.loc[df_cov[1] == row[3], 4] - row[1]) / 3

    df_cov.loc[df_cov['start'] < 0, 'start'] = 0
    df_cov['start'] = df_cov['start'].apply(lambda x: math.floor(x) if not math.isnan(x) else x)
    df_cov['end'] = df_cov['end'].apply(lambda x: math.floor(x) if not math.isnan(x) else x)

    phrase27 = get_sentence(df_cov,"UL27")
    phrase51 = get_sentence(df_cov,"UL51")
    phrase52 = get_sentence(df_cov,"UL52")
    phrase54 = get_sentence(df_cov,"UL54")
    phrase56 = get_sentence(df_cov,"UL56")
    phrase89 = get_sentence(df_cov,"UL89")
    phrase97 = get_sentence(df_cov,"UL97")
    phrase104 = get_sentence(df_cov,"UL104")
else:
    phrase27 = ""
    phrase51 = ""
    phrase52 = ""
    phrase54 = ""
    phrase56 = ""
    phrase89 = ""
    phrase97 = ""
    phrase104 = ""

# Tri des données du VCF
df = pd.read_csv(vcf_file, sep='\t', comment='#', header=None, names=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "Sample1"])
abreviations_to_letters = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}

mutations_ul = {
    27: ["A269T", "A406V", "C415S", "D534Y", "E22S", "L193F", "L335P", "L426F", "R233S", "R448P", "V353E", "W153R", "W362S"],
    51: ["P91S"],
    54: ["S290R", "D301N", "E303D", "E303G", "N408D", "N408K", "N408S", "N410K", "F412C", "F412L", "F412S", "F412V", "D413E", "D413A", "D413N", "D413Y", "K488R", "N495K", "K500N", "L501I", "T503I", "A505V", "K513E", "K513N", "K513R", "D515Y", "L516P", "L516W", "L516R", "I521T", "P522S", "P522A", "C524del", "V526L", "C539G", "C539R", "D542E", "A543P", "L545S", "L545W", "T552N", "Q578H", "Q578L", "S585A", "D588N", "F595I", "H600L", "T700A", "V715A", "V715M", "I726T", "I726V", "E756D", "E756G", "E756Q", "E756K", "L773V", "L776M", "V781I", "V787A", "V787L", "V787E", "L802M", "K805Q", "A809V", "V812L", "T813S", "T821I", "P829S", "A834P", "T838A", "G841A", "G841S", "M844T", "M844V", "V946L", "E951D", "L957F", "D981del", "L982del", "A987G"],
    56: ["C25F", "S229F", "V231A", "V231L", "N232Y", "V236A", "V236L", "V236M", "E237D", "E237G", "L241P", "T244K", "L254F", "L257I", "L257F", "K258E", "F261L", "F261C", "Y321C", "C325F", "C325Y", "C325R", "C325W", "L328V", "M329T", "A365S", "N368D", "R369M", "R369G", "R369S", "R369T"],
    89: ["N320H", "D344E", "T350M", "M359I"],
    97: ["L337M", "F342Y", "F342S", "V353A", "V356G", "K359E", "K359Q", "L397R", "L405P", "T409M", "H411Y", "H411N", "H411L", "D456N", "M460V", "M460I", "V466G", "C480F", "C480R", "C518Y", "H520Q", "P521L", "A591V", "C592G", "A594V", "A594T", "A594E", "A594G", "A594P", "A594S", "L595S", "L595F", "L595W", "E596G", "E596Y", "K599T", "C603W", "C603R", "C603S", "C607Y", "C607F", "I610T", "A613V"]
}

df['INFO_split'] = df['INFO'].str.split("|")
df['Sample1_split'] = df['Sample1'].str.split(":")

tri_type = df['INFO_split'].str[2].isin(["MODERATE", "HIGH"])
tri_freq = df['Sample1_split'].str[7].apply(lambda x: float(x.replace('%', '')) if isinstance(x, str) else x).astype(float) > 1.00
tri_ul27 = df['INFO_split'].str[3].str.contains("UL27")
tri_ul51 = df['INFO_split'].str[3].str.contains("UL51")
tri_ul52 = df['INFO_split'].str[3].str.contains("UL52")
tri_ul54 = df['INFO_split'].str[3].str.contains("UL54")
tri_ul56 = df['INFO_split'].str[3].str.contains("UL56")
tri_ul89 = df['INFO_split'].str[3].str.contains("UL89")
tri_ul97 = df['INFO_split'].str[3].str.contains("UL97")
tri_ul104 = df['INFO_split'].str[3].str.contains("UL104")

df['Acides Aminés'] = df['INFO_split'].str[10]
df['Nucléotides'] = df['INFO_split'].str[9]
df['Fréquence'] = df['Sample1_split'].str[7]
df['Profondeur'] = df['Sample1_split'].str[3]
df['taille'] = df['INFO'].str.split("|").apply(len)

df27 = df.loc[tri_type & tri_freq & tri_ul27, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df51 = df.loc[tri_type & tri_freq & tri_ul51, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df52 = df.loc[tri_type & tri_freq & tri_ul52, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df54 = df.loc[tri_type & tri_freq & tri_ul54, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df56 = df.loc[tri_type & tri_freq & tri_ul56, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df89 = df.loc[tri_type & tri_freq & tri_ul89, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df97 = df.loc[tri_type & tri_freq & tri_ul97, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
df104 = df.loc[tri_type & tri_freq & tri_ul104, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]

condition = df['taille'] >= 31
if condition.any():
      tri_type = df['INFO_split'].str[17].isin(["MODERATE", "HIGH"])
      tri_ul27 = df['INFO_split'].str[18].str.contains("UL27")
      tri_ul51 = df['INFO_split'].str[18].str.contains("UL51")
      tri_ul52 = df['INFO_split'].str[18].str.contains("UL52")
      tri_ul54 = df['INFO_split'].str[18].str.contains("UL54")
      tri_ul56 = df['INFO_split'].str[18].str.contains("UL56")
      tri_ul89 = df['INFO_split'].str[18].str.contains("UL89")
      tri_ul97 = df['INFO_split'].str[18].str.contains("UL97")
      tri_ul104 = df['INFO_split'].str[18].str.contains("UL104")
      
      df['Acides Aminés'] = df['INFO_split'].str[25]
      df['Nucléotides'] = df['INFO_split'].str[24]
      df['Fréquence'] = df['Sample1_split'].str[7]
      df['Profondeur'] = df['Sample1_split'].str[3]
      
      df27_tmp = df.loc[tri_type & tri_freq & tri_ul27, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df51_tmp = df.loc[tri_type & tri_freq & tri_ul51, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df52_tmp = df.loc[tri_type & tri_freq & tri_ul52, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df54_tmp = df.loc[tri_type & tri_freq & tri_ul54, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df56_tmp = df.loc[tri_type & tri_freq & tri_ul56, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df89_tmp = df.loc[tri_type & tri_freq & tri_ul89, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df97_tmp = df.loc[tri_type & tri_freq & tri_ul97, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df104_tmp = df.loc[tri_type & tri_freq & tri_ul104, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df27 = pd.concat([df27, df27_tmp], ignore_index=True)
      df51 = pd.concat([df51, df51_tmp], ignore_index=True)
      df52 = pd.concat([df52, df52_tmp], ignore_index=True)
      df54 = pd.concat([df54, df54_tmp], ignore_index=True)
      df56 = pd.concat([df56, df56_tmp], ignore_index=True)
      df89 = pd.concat([df89, df89_tmp], ignore_index=True)
      df97 = pd.concat([df97, df97_tmp], ignore_index=True)
      df104 = pd.concat([df104, df104_tmp], ignore_index=True)

condition = df['taille'] >= 46
if condition.any():
      tri_type = df['INFO_split'].str[32].isin(["MODERATE", "HIGH"])
      tri_ul27 = df['INFO_split'].str[33].str.contains("UL27")
      tri_ul51 = df['INFO_split'].str[33].str.contains("UL51")
      tri_ul52 = df['INFO_split'].str[33].str.contains("UL52")
      tri_ul54 = df['INFO_split'].str[33].str.contains("UL54")
      tri_ul56 = df['INFO_split'].str[33].str.contains("UL56")
      tri_ul89 = df['INFO_split'].str[33].str.contains("UL89")
      tri_ul97 = df['INFO_split'].str[33].str.contains("UL97")
      tri_ul104 = df['INFO_split'].str[33].str.contains("UL104")
      
      df['Acides Aminés'] = df['INFO_split'].str[40]
      df['Nucléotides'] = df['INFO_split'].str[39]
      df['Fréquence'] = df['Sample1_split'].str[7]
      df['Profondeur'] = df['Sample1_split'].str[3]
      
      df27_tmp = df.loc[tri_type & tri_freq & tri_ul27, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df51_tmp = df.loc[tri_type & tri_freq & tri_ul51, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df52_tmp = df.loc[tri_type & tri_freq & tri_ul52, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df54_tmp = df.loc[tri_type & tri_freq & tri_ul54, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df56_tmp = df.loc[tri_type & tri_freq & tri_ul56, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df89_tmp = df.loc[tri_type & tri_freq & tri_ul89, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df97_tmp = df.loc[tri_type & tri_freq & tri_ul97, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df104_tmp = df.loc[tri_type & tri_freq & tri_ul104, ['Acides Aminés', 'Nucléotides', 'Fréquence', 'Profondeur']]
      df27 = pd.concat([df27, df27_tmp], ignore_index=True)
      df51 = pd.concat([df51, df51_tmp], ignore_index=True)
      df52 = pd.concat([df52, df52_tmp], ignore_index=True)
      df54 = pd.concat([df54, df54_tmp], ignore_index=True)
      df56 = pd.concat([df56, df56_tmp], ignore_index=True)
      df89 = pd.concat([df89, df89_tmp], ignore_index=True)
      df97 = pd.concat([df97, df97_tmp], ignore_index=True)
      df104 = pd.concat([df104, df104_tmp], ignore_index=True)

df27["Acides Aminés"] = df27["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df51["Acides Aminés"] = df51["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df52["Acides Aminés"] = df52["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df54["Acides Aminés"] = df54["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df56["Acides Aminés"] = df56["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df89["Acides Aminés"] = df89["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df97["Acides Aminés"] = df97["Acides Aminés"].replace(abreviations_to_letters, regex=True)
df104["Acides Aminés"] = df104["Acides Aminés"].replace(abreviations_to_letters, regex=True)

df27['Acides Aminés'] = df27['Acides Aminés'].str.replace('p.', '')
df51['Acides Aminés'] = df51['Acides Aminés'].str.replace('p.', '')
df52['Acides Aminés'] = df52['Acides Aminés'].str.replace('p.', '')
df54['Acides Aminés'] = df54['Acides Aminés'].str.replace('p.', '')
df56['Acides Aminés'] = df56['Acides Aminés'].str.replace('p.', '')
df89['Acides Aminés'] = df89['Acides Aminés'].str.replace('p.', '')
df97['Acides Aminés'] = df97['Acides Aminés'].str.replace('p.', '')
df104['Acides Aminés'] = df104['Acides Aminés'].str.replace('p.', '')

df27['Nucléotides'] = df27['Nucléotides'].str.replace('c.', '')
df51['Nucléotides'] = df51['Nucléotides'].str.replace('c.', '')
df52['Nucléotides'] = df52['Nucléotides'].str.replace('c.', '')
df54['Nucléotides'] = df54['Nucléotides'].str.replace('c.', '')
df56['Nucléotides'] = df56['Nucléotides'].str.replace('c.', '')
df89['Nucléotides'] = df89['Nucléotides'].str.replace('c.', '')
df97['Nucléotides'] = df97['Nucléotides'].str.replace('c.', '')
df104['Nucléotides'] = df104['Nucléotides'].str.replace('c.', '')

if not df27.empty:
	df27["Résistance"] = df27.apply(lambda row: check_mutation_presence(row, 27), axis=1)
if not df51.empty:
	df51["Résistance"] = df51.apply(lambda row: check_mutation_presence(row, 51), axis=1)
if not df52.empty:
	df52["Résistance"] = df52.apply(lambda row: check_mutation_presence(row, 52), axis=1)
if not df54.empty:
	df54["Résistance"] = df54.apply(lambda row: check_mutation_presence(row, 54), axis=1)
if not df56.empty:
	df56["Résistance"] = df56.apply(lambda row: check_mutation_presence(row, 56), axis=1)
if not df89.empty:
	df89["Résistance"] = df89.apply(lambda row: check_mutation_presence(row, 89), axis=1)
if not df97.empty:
	df97["Résistance"] = df97.apply(lambda row: check_mutation_presence(row, 97), axis=1)
if not df104.empty:
	df104["Résistance"] = df104.apply(lambda row: check_mutation_presence(row, 104), axis=1)

# Définition d'un environnement Jinja2 avec le répertoire de templates
env = Environment(loader=FileSystemLoader(args.template_dir))
template = env.get_template('template.html')

# Transformation du DataFrame en  HTML
df27_html = df27.to_html(index=False)
df51_html = df51.to_html(index=False)
df52_html = df52.to_html(index=False)
df54_html = df54.to_html(index=False)
df56_html = df56.to_html(index=False)
df89_html = df89.to_html(index=False)
df97_html = df97.to_html(index=False)
df104_html = df104.to_html(index=False)

# Création d'un dictionnaire contexte pour le template
context = {
    'date': datetime.date.today(),
    'titre': ech,
    'sequences': sequences,
	'pourcentage': mapped_percentage,
    'tableau27': df27_html,
	'phrase27': phrase27,
    'tableau51': df51_html,
	'phrase51': phrase51,
    'tableau52': df52_html,
	'phrase52': phrase52,
    'tableau54': df54_html,
	'phrase54': phrase54,
    'tableau56': df56_html,
	'phrase56': phrase56,
    'tableau89': df89_html,
	'phrase89': phrase89,
    'tableau97': df97_html,
	'phrase97': phrase97,
    'tableau104': df104_html,
	'phrase104': phrase104
}

# Rendu du modèle Jinja2 avec les données du contexte
output = template.render(context)

# Écriture du contenu généré dans un fichier HTML
with open(output_file, 'w') as file:
    file.write(output)

print("Fichier HTML généré avec succès : output.html")
