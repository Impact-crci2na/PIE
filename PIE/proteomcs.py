import argparse
import pickle
import sys
import pyvis
from pyvis.network import Network
import networkx as nx
import matplotlib.pyplot as plt
import requests
from bioservices.uniprot import UniProt

u = UniProt(verbose = False)
import time

def convert_ACC_GENENAME(id) :
    convert = u.mapping(fr="UniProtKB_AC-ID", to="Gene_Name", query=id)
    return convert['results'][0]['to']

#Extract edges and nodes from the file
def create_initial_graph(path):
  G=nx.Graph()
  with open (path,"r") as f :
    for l in f.readlines():
        for m in l.split():
            if m not in G.nodes():
                G.add_node(m)
        G.add_edge(l.split()[0],l.split()[1])
  return G

def list_of_prot_acc(path):
  list_prot_acc=[]
  with open (path,"r") as f :
    for l in f.readlines():
        list_prot_acc.append(l[:-1])
  return list_prot_acc

def create_list_logratio(path):
  list_logratio=[]
  with open (path,"r") as f :
    for l in f.readlines():
      try:
        score=float(l.strip().replace(",","."))
        list_logratio.append(abs(score))
      except:
        print(l)
  return list_logratio

def create_list_prot_name(list_prot_acc):
    list_prot_name=[]
    dico={}
    for i in list_prot_acc:
        try:
            n=convert_ACC_GENENAME(str(i))
            list_prot_name.append(n)
            dico[n] = i
            print(n)
        except:
            print(i,"not reconised")
    return list_prot_name, dico

def create_list_proteo(list_prot_name,G):
  list_proteo=[]
  for i in list_prot_name:
    if i in G.nodes():
        list_proteo.append(i)
    else:
      print(i,"not in the databases")
  return(list_proteo)



def create_dico_logratio(list_prot_name,list_logratio):
  dico_logratio={}
  for i in range(len(list_logratio)):
    try:
      dico_logratio[list_prot_name[i]]=list_logratio[i]
    except:
      print(i)
  return(dico_logratio)


def display_shortest_path_static(start_GN, list_prot, G, dico_logratio):
    print("Start pathfindddding...")

    # Charger le graphe à partir d'un fichier pickle

    # Initialisation d'un nouveau graphe NetworkX pour les chemins les plus courts
    G1 = nx.Graph()

    # Ajouter les chemins les plus courts au graphe G1
    for prot in list_prot:
        try:
            path = nx.shortest_path(G, source=start_GN, target=prot)
            nx.add_path(G1, path)
        except nx.NetworkXNoPath:
            print(f"{prot} does not have interactions with {start_GN}")

    # Préparation de la figure
    # plt.figure(figsize=(100, 40))
    # pos = nx.spring_layout(G1)  # positions for all nodes
    #
    # # Nodes drawing
    # nx.draw_networkx_nodes(G1, pos, node_color='green', node_size=50)
    # nx.draw_networkx_nodes(G1, pos, nodelist=[start_GN], node_color='red', node_size=100)
    # nx.draw_networkx_nodes(G1, pos, nodelist=list_prot, node_color='blue', node_size=[dico_logratio.get(node, 1) * 100 for node in list_prot])
    #
    # # Edges drawing
    # nx.draw_networkx_edges(G1, pos, width=1.0, alpha=0.5)
    #
    # # Labels drawing
    # nx.draw_networkx_labels(G1, pos, font_size=10, font_color='black')
    color_map=[]
    for node in G1:
        if node==start_GN:
            color_map.append('red')
        elif node in list_prot:
            color_map.append("blue")
        else:
            color_map.append('green')
    size=[]
    for n in G1.nodes():
        if n==start_GN:
            size.append(50000)
        elif n in dico_logratio.keys():
            size.append(dico_logratio[n]*3000)
        else:
            size.append(1000)
    fig = plt.figure(1, figsize=(100, 40))
    pos = nx.spring_layout(G1,scale=4)
    nx.draw(G1,pos,with_labels=True,node_color=color_map,node_size=size,bbox=dict(facecolor="white", edgecolor='black', boxstyle='round,pad=0.2'), font_size=25)
    # Sauvegarder le graphique dans un fichier image
    plt.savefig('mygraph.png')
    plt.close()

    return G1


def display_shortest_path(start_GN, list_prot, G, dico_logratio):
    print("Start pathfinding...")

    # Initialisation d'un nouveau graphe NetworkX pour les chemins les plus courts
    G1 = nx.Graph()

    # Ajouter les chemins les plus courts au graphe G1
    for prot in list_prot:
        try:
            path = nx.shortest_path(G, source=start_GN, target=prot)
            nx.add_path(G1, path)
        except nx.NetworkXNoPath:
            print(f"{prot} don't have interactions with {start_GN}")

    # Création du réseau Pyvis avec configuration des options de physique
    net = Network(notebook=False, height="750px", width="100%", bgcolor="#222222", font_color="white")
    net.barnes_hut(gravity=-8000, central_gravity=0.3, spring_length=100, spring_strength=0.01, damping=0.09, overlap=0)

    # Configurer les attributs des nœuds directement dans Pyvis
    for node in G1.nodes():
        color = 'red' if node == start_GN else 'blue' if node in list_prot else 'green'
        default_size = 3 if node not in list_prot else dico_logratio.get(node, 1) * 5
        size = 15 if node == start_GN else default_size
        title = node  # Titre pour affichage au survol
        net.add_node(node, label=node, color=color, size=size, title=title)

    # Importer les chemins (arêtes) depuis G1
    for edge in G1.edges():
        net.add_edge(*edge)

    # Affichage du graphique
    net.show('mygraph.html')
    return G1, net

def top_node(G1):
  top={}

  for node, degree in G1.degree():
      if degree not in top:
          top[degree] = [node]
      else:
          top[degree].append(node)
  return top

def top_connexion(G1):
    connections = {}
    for node in G1.nodes():
        # Récupère tous les voisins (noeuds connectés) du noeud actuel
        neighbors = list(G1.neighbors(node))
        connections[node] = neighbors
    return connections

def print_top_connexion(dico_degree, top):
    # Trier les clés par ordre décroissant
    classement = sorted(dico_degree.keys(), reverse=True)
    print(f"top {top} : ")
    # S'assurer de ne pas dépasser le nombre de clés disponibles
    top_effectif = min(top, len(classement))
    for i in range(top_effectif):
        print(f"{classement[i]} : {dico_degree[classement[i]]}")


def save_graph(graph, filename):
    with open(filename, 'wb') as f:
        pickle.dump(graph, f)

def fetch_uniprot_ids_with_bioservices(gene_name, organism="Homo sapiens"):
    """
    Fetch UniProt IDs directly using the UniProt API.

    :param gene_name: str, the gene name to query for.
    :param organism: str, the organism (species) of the gene.
    :return: list of UniProt IDs
    """
    url = f"https://rest.uniprot.org/uniprotkb/search?query=organism_id:9606+AND+gene_exact:\"{gene_name}\"&format=tsv"
    response = requests.get(url)
    if response.status_code == 200:
        lines = response.text.strip().split('\n')
        if len(lines) > 1:  # Vérifier s'il y a au moins une ligne de données
            first_data_line = lines[1]  # La ligne après les en-têtes
            first_id = first_data_line.split('\t')[0]  # Récupérer l'identifiant de la première colonne
            return first_id
    else:
        print("Failed to retrieve data: HTTP", response.status_code)
        return []


def write_node_on_txt(G,filename="complete_connexion"):
    with open(f"{filename}.txt", 'w') as fichier:
        for noeud in G.nodes():
            fichier.write(f"{noeud}\n")
def load_graph(file_path):
    """
    Charge un graphe NetworkX à partir d'un fichier pickle.
    """
    try:
        with open(file_path, 'rb') as file:
            return pickle.load(file)
    except EOFError:
        print(f"Warning: Empty or corrupted file encountered: {file_path}")
        return None
def display_shortest_path_staticx(start_GN, list_prot, G1, dico_logratio):
    print("Start pathfinding...")
    color_map=[]
    for node in G1:
        if node==start_GN:
            color_map.append('red')
        elif node in list_prot:
            color_map.append("blue")
        else:
            color_map.append('green')
    size=[]
    for n in G1.nodes():
        if n==start_GN:
            size.append(50000)
        elif n in dico_logratio.keys():
            size.append(dico_logratio[n]*3000)
        else:
            size.append(1000)
    fig = plt.figure(1, figsize=(60, 40))
    pos = nx.spring_layout(G1,scale=4)
    nx.draw(G1,pos,with_labels=True,node_color=color_map,node_size=size,bbox=dict(facecolor="white", edgecolor='black', boxstyle='round,pad=0.2'), font_size=25)
    # Sauvegarder le graphique dans un fichier image
    plt.savefig('mygraph.png')
    plt.close()
def load_graph(file_path):
    """
    Charge un graphe NetworkX à partir d'un fichier pickle.
    """
    try:
        with open(file_path, 'rb') as file:
            return pickle.load(file)
    except EOFError:
        print(f"Warning: Empty or corrupted file encountered: {file_path}")
        return None



def list_of_intermediate(graph, list_tot):
    network = set(graph.nodes)
    protein_inter = network.difference(list_tot)
    with open("intermediate.txt", "w") as file:
        for proteins in protein_inter:
            file.write(proteins + "\n")
