# -*- coding: utf-8 -*-
"""s30338_zadanie1.ipynb

Zadanie 1: Wizualizacja danych ekspresji genów

Wykorzystaj biblioteki pandas, seaborn i matplotlib do analizy danych ekspresji genów z badania porównującego komórki nabłonka dróg oddechowych u palaczy i osób niepalących
(dataset GDS2490.soft.gz : https://drive.google.com/file/d/1tb_imtWMufVWXuOJtiXR6eFwPD6Wlp1T/view?usp=sharing ).


Dane znajdują się w pliku GDS2490.soft i zawierają wartości ekspresji tysięcy genów w 11 próbkach:

    GSM114084-GSM114088: non-smoker (kontrola, n=5)
    GSM114078-GSM114083: smoker (palący, n=6)

Napisz skrypt w Pythonie z następującymi wymaganiami:
1. Wczytaj dane z pliku GDS2490.soft
    Pomiń linie zaczynające się od # lub ^
    Wczytaj tabelę zaczynającą się od ‘!dataset_table_begin’ i kończy na: '!dataset_table_end'
    Użyj kolumny IDENTIFIER jako nazwy genów
2. Przekształć dane do formatu "long" (funkcja pandas.melt)
    Kolumny: gene, sample, expression, group
    group = "non-smoker" lub "smoker" (na podstawie ID próbki)
3. Zidentyfikuj 10 genów o największej różnicy ekspresji między grupami
    Oblicz średnią ekspresję dla każdej grupy
    Oblicz różnicę absolutną |mean_smoker - mean_non_smoker|
4. Wybierz top 10 genów
    Dla tych 10 genów wyświetl podstawowe statystyki:
    
    - Nazwy genów
    - Średnia ekspresja w grupie non-smoker
    - Średnia ekspresja w grupie smoker
    - Różnica (fold change)
5. Utwórz trzy wykresy (tylko dla wybranych 10 genów).
    - boxplot – porównanie ekspresji między grupami dla każdego genu      (x=gene, y=expression, hue=group)
    - violinplot – rozkład z gęstością, split by group
    - heatmap – ekspresja 10 genów w 11 próbkach (wiersze=geny, kolumny=próbki pogrupowane: non-smoker | smoker)

Możesz zainspirować się formatami z galerii seaborn i wybrać inne wizualizacje.

##============================================================================

Przykładowy output:

=== ANALIZA EKSPRESJI GENÓW ===

Wczytano dane: 1000 genów, 11 próbek

Grupy: non-smoker (5), smoker (6)

Top 10 genów z największą różnicą ekspresji:
    
    Gene        Non-smoker  Smoker   Difference
    CYP1B1      150.2       1250.8   +1100.6
    ALDH3A1     89.5        650.3    +560.8
    NQO1        120.4       550.1    +429.7
    ...

Zapisano wykresy:
    
    - boxplot_top10_genes.png
    - violinplot_top10_genes.png
    - heatmap_top10_genes.png

Przygotuj finalną wersję skryptu z kodem Python (plik .py), który będzie rozwiązaniem dla tego zadania.
"""

#todo dopracowac output na py

#importy

#pip install pandas seaborn matplotbil numpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

from io import StringIO

"""Wczytaj dane z pliku GDS2490.soft
    
    Pomiń linie zaczynające się od # lub ^
    Wczytaj tabelę zaczynającą się od ‘!dataset_table_begin’ i kończy na: '!dataset_table_end'
    Użyj kolumny IDENTIFIER jako nazwy genów
"""
print("=== ANALIZA EKSPRESJI GENÓW ===")

sciezka_do_pliku = "GDS2490.soft"

table_lines = []
inside_table = False

with open(sciezka_do_pliku, 'r') as file:
  for line in file:
    line = line.strip()

    #poczatek tabeli
    if line == "!dataset_table_begin":
      inside_table = True
      continue

    #koniec tabeli
    if line == "!dataset_table_end":
      inside_table = False
      break

    #w tabeli
    if inside_table:
      #pominiecie linii zaczynających się od # lub ^
      if line.startswith("#") or line.startswith("^"):
        continue

      table_lines.append(line)

# print(table_lines)

#nazwy genow:
gen_names = []
for gen in table_lines[0].split("\t")[2:]:
  gen_names.append(gen)

"""2. Przekształć dane do formatu "long" (funkcja pandas.melt)

    Kolumny: gene, sample, expression, group
    
    group = "non-smoker" lub "smoker" (na podstawie ID próbki)
"""

rows = [linie.split("\t") for linie in table_lines]

#naglowki
header = rows[0]
dane = rows[1:]

df = pd.DataFrame(dane, columns=header)
df_long = pd.melt(
    df,
    id_vars=["ID_REF","IDENTIFIER"],
    var_name="sample",
    value_name="expression"
)

#zmiana nazwy Idendifier -> gene
df_long = df_long.rename(columns={"IDENTIFIER":"gene"})

non_smokers = [
    "GSM114084", "GSM114085", "GSM114086",
    "GSM114087", "GSM114088"
]

smokers = [
    "GSM114078", "GSM114079", "GSM114080",
    "GSM114081", "GSM114082", "GSM114083"
]

df_long["group"] = df_long["sample"].apply(
    lambda x: "smoker" if x in smokers else "non-smoker"
)

print(f"Wczytano dane: {len(df_long["gene"])} genów, {len(df_long["sample"])} próbek")
print(f"Grupy: non-smoker {len(non_smokers)}, smoker {len(smokers)}\n")

# expression -> float
df_long["expression"] = pd.to_numeric(df_long["expression"])

"""3. Zidentyfikuj 10 genów o największej różnicy ekspresji między grupami
    
    Oblicz średnią ekspresję dla każdej grupy
    
    Oblicz różnicę absolutną
      |mean_smoker - mean_non_smoker|
"""

#srednia ekspresja dla kazdego genu
mean_expression = (
    df_long
    .groupby(["gene","group"])["expression"]
    .mean()
    .unstack()
)

#roznica absolutna
mean_expression["abs_difference"] = abs(
    mean_expression["smoker"] - mean_expression["non-smoker"]
)

#fold
mean_expression["fold_difference"] = (
    mean_expression["smoker"] / mean_expression["non-smoker"]
)

#10 genow o najwiekszej roznicy
top10genes = (
    mean_expression
    .sort_values(by="abs_difference", ascending=False)
    .head(10)
)
"""4. Wybierz top 10 genów
    Dla tych 10 genów wyświetl podstawowe statystyki:
    
    - Nazwy genów
    - Średnia ekspresja w grupie non-smoker
    - Średnia ekspresja w grupie smoker
    - Różnica (fold change)
"""

top10genes_stats = top10genes[
    ["smoker","non-smoker","fold_difference","abs_difference"]
]

print("Top 10 genów z największą różnicą ekspresji:",top10genes_stats, sep="\n")

"""5. Utwórz trzy wykresy (tylko dla wybranych 10 genów).
    - boxplot – porównanie ekspresji między grupami dla każdego genu      (x=gene, y=expression, hue=group)
    - violinplot – rozkład z gęstością, split by group
    - heatmap – ekspresja 10 genów w 11 próbkach (wiersze=geny, kolumny=próbki pogrupowane: non-smoker | smoker)
"""

top10genes_names = top10genes.index.tolist()
df_top10genes = df_long[df_long["gene"].isin(top10genes_names)]

print("Zapisano wykresy:")

#boxplot
plt.figure(figsize=(10,6))

sns.boxplot(
    data=df_top10genes,
    x="gene",
    y="expression",
    hue="group"
)

plt.xticks(rotation = 45)
plt.title("Expression boxplot for top 10 genes")
plt.xlabel("Gene")
plt.ylabel("Expression")
plt.tight_layout()
plt.savefig("boxplot_top10_genes.png")
plt.show()
print("- boxplot_top10_genes.png")

#violinplot
plt.figure(figsize=(10,6))
sns.violinplot(
    data=df_top10genes,
    x="gene",
    y="expression",
    hue="group",
    split=True
)

plt.xticks(rotation = 45)
plt.title("Expression violinplot for top 10 genes")
plt.xlabel("Gene")
plt.ylabel("Expression")
plt.tight_layout()
plt.savefig("violinplot_top10_genes.png")
plt.show()
print("- violinplot_top10_genes.png")

#heatmap
heatmapData = (
    df_top10genes
    .pivot_table(
        index="gene",
        columns="sample",
        values="expression"
    )
)

orderedCollumns = non_smokers + smokers
heatmapData = heatmapData[orderedCollumns]

plt.figure(figsize=(10,6))
sns.heatmap(
    heatmapData,
    cmap= "viridis"
)

plt.title("Expression heatmap for top 10 genes")
plt.xlabel("Sample")
plt.ylabel("Gene")
plt.tight_layout()
plt.savefig("heatmap_top10_genes.png")
plt.show()
print("- heatmap_top10_genes.png")