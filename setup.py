from setuptools import setup, find_packages

with open('README.md') as f:
    long_description = f.read()

setup(name='txfeature',
      version='1.0',
      packages=find_packages(),
      install_requires=['pybedtools', 'Bio', 'pandas', 'biopython'],

      # metadata to display on PyPI
      author="Waqar Arif",
      author_email="me@example.com",
      description="A tool to aggregate transcript features from gene annotation files.",
      license='MIT',
      long_description=long_description,
      keywords="bioinformatics rna-seq gene expression",
      url="https://github.com/warif2/txfeature",
      project_urls={
            "Bug Tracker": "Coming Soon",
            "Documentation": "Coming Soon"
      }
      )
