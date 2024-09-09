from setuptools import setup, find_packages

setup(
    name="Proteomics Interactions Expressions",  # Remplace par le nom de ton package
    version="1.1",
    packages=find_packages(),  # Recherche automatiquement les packages dans ton répertoire
    install_requires=[
        "argparse",      # Pour l'analyse des arguments de ligne de commande
        "pickle-mixin",  # Pour la gestion des objets Pickle
        "pyvis",         # Pour la visualisation des graphes
        "networkx",      # Pour la manipulation des graphes
        "matplotlib",    # Pour les graphiques et visualisations
        "requests",      # Pour faire des requêtes HTTP
        "bioservices"    # Pour accéder aux services web comme UniProt
    ],
    description="A package for network visualization and bioinformatics",
    author="Karen Sobriel, Grégoire Menard, Haladi Ayad",
    author_email="crcina.impact.bioinfo@gmail.com",
    url="https://gitlab.univ-nantes.fr/E179974Z/pie.git",  # Lien vers ton dépôt GitHub
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',  # Version de Python minimum requise
)
