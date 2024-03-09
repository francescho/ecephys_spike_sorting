
import pandas as pd
import sys

def main():
	# path to session attributes: read from 2nd argument from command line 
	_, attributes_path = sys.argv

	# read in session attributes
	sbatch_attributes = pd.read_csv(attributes_path)
	print(sbatch_attributes)
	
if __name__ == "__main__":
    main()
