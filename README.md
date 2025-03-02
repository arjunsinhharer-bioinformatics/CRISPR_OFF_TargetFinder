To run the script download the zip file, open the file, and then save the file to your desktop.
The name of the file is “Arjunsinh_Harer_Final_Project_RBIF_109”.
Navigate to this folder in your terminal and type the following command.
python o _target.py
You will now see an example output of the script. The script evaluates the guide suitability
of guideRNAs needed to knock out EFGR in homo-sapiens.
This script allows you to evaluate guideRNA candidates for a CRISPR knockout experiment
in homo-sapiens.
If you want to run the script on a particular gene of interest, you can start by running the
following command.
python download_fasta.py “Gene Name of Choice” “Your email”
An example command looks like this:
python download_fasta.py EGFR arjunsinhharer@brandeis.edu
Open the “Arjunsinh_Harer_Final_Project_RBIF_109” directory and ensure your FASTA for
your gene of interest is downloaded.
After downloading the file delete the any existing txt files in the directory if they exist. The
o _target.py script can only run if there is one .txt file in the directory.
This is the only manual step of the process. You need to obtain a list of guide RNAs for your
gene of interest. You can do this using the Synthego Knockout Guide Design Tool
Once you open the link you’ll see a web-based tool to generate guide RNAs for a gene. For
genome type Homo sapiens, select any gene of your choice and then click “search”.
Copy paste the recommended guide sequences in a .txt file, and each line should have one
guide sequence listed. Once you copy-paste your guide sequences into the text file, save it
to the directory “Arjunsinh_Harer_Final_Project_RBIF_109”.
Please name the text file with the following convention otherwise the script won’t run
guides_ {Gene name of your choice}.txt.
Now that you have the FASTA file for your gene of interest as well as the guide RNA
sequences you are now ready to run the script and get the output.
You can do this by simply typing the command specified in the start of this readme file
which is.
python o _target.py
Great! You have now successfully run your own independent analysis of guide RNA
suitability for a CRISPR experiment. If you run into any issues running this script, please do
not deduct points, instead call me, or text me 703-582-3705. You can text this number as
well, or email me at arjunsinhharer@brandeis.edu Thanks, and I hope you enjoy!
