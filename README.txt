CoVaMa Version 0.1
Last Modified: 4-8-15

Contents:

	TestData_FHV-25k (15 files)
		Small test data set using 25000 reads from an FHV RNAseq experiment
		- FHV_Genome_padded.?.ewbt (6 files, bowtie index)
		- FHV_R1_25k.txt	(Raw Reads_5'read)
		- FHV_R2_25k.txt	(Raw Reads_3'read)
		- FHV_R1_25k_FHV-mapping.sam	(FHV Genome mapping, 5' reads)
		- FHV_R2_25k_FHV-mapping.sam	(FHV Genome mapping, 3' reads)
		- R1_ViReMa_Output.txt	(No Header)
		- R2_ViReMa_Output.txt	(No Header)
		- R1_mFHV_mDm.txt	(5' reads unaligned to FHV and D.melanogaster genome)
		- R1_mFHV_mDm.txt	(3' reads unaligned to FHV and D.melanogaster genome)
		- FHV_Genome_padded.txt	Contains reference genes for Flock House Virus with long 3' terminal A residues.


	CoVaMa_Make_Matrices.py
		Script used to generate Nucleotide Matrices.  Runs from Command-line.

	CoVaMa_Make_Matrices_wAA.py
		Script used to generate Amino Acid Matrices.  Runs from Command-line.
	
	CoVaMa_Merge_Matrices.py
		Script used to Merge multiple matrices.  Runs from Command-line.

	CoVaMa_Analyse_Matrices.py
		Script used to Analyse Nucleotide  Matrices.  Runs from Command-line.
	
	CoVaMa_Analyse_Matrices_wAA.py
		Script used to Analyse Amino Acid Matrices.  Runs from Command-line.

	Config.py
		This script carries the global variables.

	README.txt
		Includes instructions to run CoVaMa.

	LICENSE.txt
		Copyright info.


CoVaMa is a simple python script and so should not require any special installation.  

CoVaMa requires python version 2.7 and Numpy. 

To analyse recombination/fusion events, CoVaMa requires the output files from ViReMa_0.7

CoVaMa is run from the command line:

>python CoVaMa_Make_Matrices.py Virus_Index Input_Data Output_Directory [args]
>python CoVaMa_Analyse_Matrices.py PickleFile Output_Data Output_Directory [args]

If amino acid sequences wish to be analyse rather than nucleotide sequences, use the '_wAA' versions of the scripts.
Open reading frames must be given in the command-line for CoVaMa_Make_Matrices_wAA.py
Currently, only one open reading frame can be used, and only one gene can be analysed. 
Split the alignment data accrodingly if multiple genes are to be analysed.


Example using test data:

>python CoVaMa_Make_Matrices.py TestData TestData/FHV_Genome_padded.txt --SAM1 TestData/FHV_R1_25k_FHV-mapping.sam --ViReMa_Output TestData/R1_ViReMa_Output.txt --Min_Fusion_Coverage 10 --SAM2 TestData/FHV_R2_25k_FHV-mapping.sam --ViReMa_Output2 TestData/R2_ViReMa_Output.txt

>python CoVaMa_Analyse_Matrices.py TestData/Total_Matrices.py.pi Output.txt TestData/  --Min_Coverage 10 --Min_Fusion_Coverage 10 -OutArray -Weighted

------------------------------------------------------------------------------

CoVaMa_Make_Matrices.py 
			[-h] [--SAM1 SAM1] [--SAM2 SAM2]
                        [--ViReMa_Output VIREMA_OUTPUT]
                        [--ViReMa_Output2 VIREMA_OUTPUT2]
                        [--Min_Coverage_Output MIN_COVERAGE_OUTPUT]
                        [--Fusion_Exclusion FUSION_EXCLUSION]
                        [--Min_Fusion_Coverage MIN_FUSION_COVERAGE]
                        [--Ends ENDS]
                        Output_Dir FASTAFile

Required arguments:

  Output_Dir            Enter name of desired output directory

  FASTAFile             Enter name of Fasta reference genome used in original
                        sequence alignment e.g. FHV_Genome.txt

Optional arguments:

  --SAM1 SAM1           Enter name of SamFile

  --SAM2 SAM2           Enter name of second/paired SamFile

  --ViReMa_Output VIREMA_OUTPUT
                        Enter name of second/paired SamFile

  --ViReMa_Output2 VIREMA_OUTPUT2
                        Enter name of second/paired SamFile

  --Fusion_Exclusion FUSION_EXCLUSION
                        Required number of nucleotides to exclude
                        recombination event. Default == 10

  --Min_Fusion_Coverage MIN_FUSION_COVERAGE
                        Enter minimum coverage over pairs of associated
                        nucleotides. Default value is 1000

  --Ends ENDS           Enter number of nucleotides to ignore from 5' and 3'
                        extremities. Default value is 0


------------------------------------------------------------------------------

CoVaMa_Make_Matrices_wAA.py 
			[-h] [--SAM1 SAM1] [--SAM2 SAM2]
                        [--ViReMa_Output VIREMA_OUTPUT]
                        [--ViReMa_Output2 VIREMA_OUTPUT2]
                        [--Min_Coverage_Output MIN_COVERAGE_OUTPUT]
                        [--Fusion_Exclusion FUSION_EXCLUSION]
                        [--Min_Fusion_Coverage MIN_FUSION_COVERAGE]
                        [--Ends ENDS]
                        Output_Dir FASTAFile

Required arguments:

  Output_Dir            Enter name of desired output directory

  FASTAFile             Enter name of Fasta reference genome used in original
                        sequence alignment e.g. FHV_Genome.txt

  ORFStartNuc           Enter minimum coverage over pairs of associated
                        nucleotides. Default value is 1

  ORFFinishNuc          Enter minimum coverage over pairs of associated
                        nucleotides. Default value is 10000


Optional arguments:

  --SAM1 SAM1           Enter name of SamFile

  --SAM2 SAM2           Enter name of second/paired SamFile

  --ViReMa_Output VIREMA_OUTPUT
                        Enter name of second/paired SamFile

  --ViReMa_Output2 VIREMA_OUTPUT2
                        Enter name of second/paired SamFile

  --Fusion_Exclusion FUSION_EXCLUSION
                        Required number of nucleotides to exclude
                        recombination event. Default == 10

  --Min_Fusion_Coverage MIN_FUSION_COVERAGE
                        Enter minimum coverage over pairs of associated
                        nucleotides. Default value is 1000

  --Ends ENDS           Enter number of nucleotides to ignore from 5' and 3'
                        extremities. Default value is 0

--------------------------------------------------------------------------

CoVaMa_Merge_Matrices.py 
			[-h] [--Multiplier MULTIPLIER]
                        [--Min_Coverage MIN_COVERAGE]
                        PickleFiles Output_File Output_Dir

Required arguments:

  PickleFiles           Enter name of PickleFiles as string 			(e.g.: 'Total_Matrices1.py.pi Total_Matrices2.py.pi')

  Output_File           Enter name of desired output file

  Output_Dir            Enter name of desired output directory

Optional arguments:

  --Multiplier MULTIPLIER
                        Enter integer number of reads desired per table,
                        Default = 1'000'000

  --Min_Coverage MIN_COVERAGE
                        Number of mapped nucleotides per contigency table
                        required for merge. Default == 100

---------------------------------------------------------------------------

CoVaMa_Analyse_Matrices.py 
			[-h] [--Min_Coverage MIN_COVERAGE]
                        [--Min_Fusion_Coverage MIN_FUSION_COVERAGE]
                        [--Min_Pop_Fraction MIN_POP_FRACTION]
                        [-OutArray] [-OutAllLDs] [-Merge]
                        [-Weighted]
                        [--Min_Merged_Matrices MIN_MERGED_MATRICES]
                        [--LD_Precision LD_PRECISION]
                        PickleFile Output_File Output_Dir

Required arguments:

  PickleFile            Enter name of PickleFile

  Output_File           Enter name of desired output file

  Output_Dir            Enter name of desired output directory

Optional arguments:

  --Min_Coverage MIN_COVERAGE
                        Enter minimum coverage over pairs of associated
                        nucleotides. Default value is 1000

  --Min_Fusion_Coverage MIN_FUSION_COVERAGE
                        Enter minimum coverage over recombination event.
                        Default value is 1000

  --Min_Pop_Fraction MIN_POP_FRACTION
                        Enter a float between 0 and 1 for the required
                        association fraction. Default value is 0.001.

  -OutArray             Write array values into output file if LD is found.
                        Default value is false.

  -OutAllLDs            Write out all LD values found. Default value is output
                        highest value only.

  -Merge                If input Matrix is from merged matrices, select
                        -Merge. Default is False

  -Weighted             Weight each LD value according to fraction of reads
                        from total contingency table. Default is False.

  --Min_Merged_Matrices MIN_MERGED_MATRICES
                        Enter a minimum number of matrices required that were
                        used to generated merged matrix. Default value is 1.

  --LD_Precision LD_PRECISION
                        Enter number of floating point digits for LD
                        output. Default value is 6.

-------------------------------------------------------------------------------

CoVaMa_Analyse_Matrices_wAA.py 
			[-h] [--Min_Coverage MIN_COVERAGE]
                        [--Min_Fusion_Coverage MIN_FUSION_COVERAGE]
                        [--Min_Pop_Fraction MIN_POP_FRACTION]
                        [-OutArray] [-OutAllLDs] [-Merge]
                        [-Weighted]
                        [--Min_Merged_Matrices MIN_MERGED_MATRICES]
                        [--LD_Precision LD_PRECISION]
                        PickleFile Output_File Output_Dir

Required arguments:

  PickleFile            Enter name of PickleFile

  Output_File           Enter name of desired output file

  Output_Dir            Enter name of desired output directory

Optional arguments:

  --Min_Coverage MIN_COVERAGE
                        Enter minimum coverage over pairs of associated
                        nucleotides. Default value is 1000

  --Min_Fusion_Coverage MIN_FUSION_COVERAGE
                        Enter minimum coverage over recombination event.
                        Default value is 1000

  --Min_Pop_Fraction MIN_POP_FRACTION
                        Enter a float between 0 and 1 for the required
                        association fraction. Default value is 0.001.

  -OutArray             Write array values into output file if LD is found.
                        Default value is false.

  -OutAllLDs            Write out all LD values found. Default value is output
                        highest value only.

  -Merge                If input Matrix is from merged matrices, select
                        -Merge. Default is False

  -Weighted             Weight each LD value according to fraction of reads
                        from total contingency table. Default is False.

  --Min_Merged_Matrices MIN_MERGED_MATRICES
                        Enter a minimum number of matrices required that were
                        used to generated merged matrix. Default value is 1.

  --LD_Precision LD_PRECISION
                        Enter number of floating point digits for LD
                        output. Default value is 6.

-------------------------------------------------------------------------------