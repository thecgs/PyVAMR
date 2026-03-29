from setuptools import setup, find_packages

setup(name="pyvamr",
      version="1.0.0",
      description='A Python package for visualizing animal mitochondrial rearrangements',
      url='https://github.com/thecgs/PyVAMR',
      author='Guisen Chen',
      author_email='thecgs001@foxmail.com',
      long_description="A Python package for visualizing animal mitochondrial rearrangements",
      license='MIT License',
      packages=find_packages(exclude=["doc", ".github"]),
      keywords=['bioinformatics', 'mitochondrial', 'genome', 'rearrangement'],
      install_requires=["biopython",
                        "matplotlib"]
     )
