use warnings;
use strict;
use List::Util qw(min max); # For Finding Max and Min in a list

# Syntax for using this program :
# perl AnalyzeLipidomes.pl <FolderName>
# Example : perl AnalyzeLipidomes.pl Test

my $MethodChoice = 4;

# By default, Levenshtein Distance is used to determine similarity betweeen molecules. 
# It is possible to use other methods also. To do so, syntax must be changed to..

# perl 121012_isomers.pl <FolderName> <MethodNumber>
# Example : perl 121012_isomers.pl 121012_YeastTest 4

# If you are changing method, un-comment next line
# $MethodChoice = $ARGV[1];

# MethodNumber Choices :
# 1. Fingerprint method
# 2. Word Frequency (LINGO) method
# 3. Bioisosteric method
# 4. Levenshtein distance
# 5. (Smith-Waterman) Local Alignment
# 6. (MUSCLE) Multiple Sequence Alignment

# Saving Input Folder Name
my $InputFolderName = $ARGV[0];

# Along with SMILES containing files, input folder must contain two other files :
# 1. with name "UserSMILES" : A plain text file that must contain SMILES for those lipid species that are not supported by LipidMapsTools. If you are not sure which Lipid species are not supported by LipidMapsTools, run this program with empty "UserSMILES" file and this will program will tell you which lipid species are not supported by LipidMapsTools
# 2. with name "FA" : A plain text file that will tell the program which fatty acids to use.
# These two files must follow the same formatting as that of example files.

# Changes from previous version of this program : 
# 1. Will take multiple files as input (All files must be placed in a folder and folder name is taken as input from user)
# 2. Will read concentrations from input files

# Caveats :
# 1. This program uses bash commands, so this program works only on linux OS
# 2. The format of molecule abbreviations in input file must be exactly as they are in example files. For instance, PI 18:1-18:1 has to be written as 
#    PI 18:1/18:1. Similarly other lipid species also have certain fixed abbreviation formats. This program will read a molecule only if it is given
#    in correct format !
# 3. This program is meant to generate structure possibilities for lipid species that have defined SN1 and SN2 Chains (Eg. PI 18:1/16:2)
#    For instance, this program cannot generate isomer possibilities for PI 36:2. 
#    If you donot have SN1 and SN2 information, use 120418_isomers.pl instead !
#    /media/second HDD 1/chakravarthy_backup/DesktopDump120717/Dm_lipidome/120418_isomers.pl

# Declaring Variables
my $VarString1;
my $VarString2;
my @VarList1;
my @VarList2;
my ($c1,$c2,$c3,$c4,$d1,$d2,$d3,$d4,$h1,$h2,$h_position,$db_positions1,$db_positions2,$db_positions3,$db_positions4,$i);

my @ClassList=(); # Initializing an array to store lipid classes
my @CentroidList=();
my @FAVarList2=(); # Probably not required. Included for Debugging
my $LipidName;
my $LipidNameCopy;
my $centroid;
my $SMILESlistFileName;
my $CentroidListFileName;
my @LipidNamesList=();
my $LipidClass;
my $LipidMapsProgramToUse;
my @Chains; # A Temporary list to store fatty acid chains
my @UnsupportedSpecies; # Array to store unused lipid abbreviations 

my $LysoFlag=0; # If 0, its NOT a lyso species. If 1, its a lyso species.

# Obtaining Names of files present in InputFolder
my @FileNameList = `ls ./$InputFolderName`;
if (scalar(@FileNameList) == 0){ die "\nSomething is wrong, didnt find any files in $InputFolderName folder. Program terminating\n\n";}

# Removing known FileNames ("UserSMILES" and "FA") from @FileNameList
my $FileNameList = join(" ", @FileNameList);
$FileNameList =~ s/UserSMILES//g;
$FileNameList =~ s/FA//g;
@FileNameList = split (/\s+/, $FileNameList);

# Checking for duplicate FileNames in InputFolder
my @UniqFileNameList = uniq(@FileNameList); 
if (scalar(@UniqFileNameList) != scalar(@FileNameList)) { 
  die "\nThere are ",  scalar(@FileNameList) - scalar(@UniqFileNameList), " duplicate file names in your input folder.... program is terminating\n\n";
}

# Saving number of Input files
my $NumOfInputFiles = scalar(@FileNameList);

# Reading UserSMILES file
my %UserSMILES; # Hash table to store UserSMILES
open UserSMILES_FILE, "./$InputFolderName/UserSMILES";
while (<UserSMILES_FILE>){
@_ = split (/\t/,$_);
chomp($_[1]);
$UserSMILES{$_[1]}="$_[0]";
}
close UserSMILES_FILE;

# Reading each file in InputFolder and making a combined list of lipid species
my $TempFileName;
my $HighestConc=0; # Variable to store highest concentration observed in given set of files. Useful for plotting.
my @LipidsList; # A list to store all lipid species
my @TempLipidsVarList1;
my @TempLipidsVarList2;

foreach $TempFileName (@FileNameList){

	@TempLipidsVarList1=();
	@TempLipidsVarList2=();
	open InFile, "./$InputFolderName/$TempFileName" or die "Can't open file $TempFileName in folder $InputFolderName : $!";
	while (<InFile>){
		@_ = split(/\t/,$_); # Separating LipidSpecies Names from their Concentrations
		push @TempLipidsVarList1, $_[0];
		chomp($_[1]);
		if ($HighestConc < $_[1]) { $HighestConc = $_[1];}
	}
	close InFile;

	# Checking if there are any duplicate entries in a given file.
	@TempLipidsVarList2 = uniq(@TempLipidsVarList1);
	if (scalar(@TempLipidsVarList2) != scalar(@TempLipidsVarList1)) { 
		die "\nThere are ",  scalar(@TempLipidsVarList1) - scalar(@TempLipidsVarList2), " duplicate molecule names in file $TempFileName. Program terminating\n\n";
	}
	
	push @LipidsList,@TempLipidsVarList1;
}

# Combining lipid species from all files into one single Duplicate-Free list.
@TempLipidsVarList2 = ();
@TempLipidsVarList2 = uniq(@LipidsList);
@LipidsList = ();
@LipidsList = @TempLipidsVarList2;

# Giving path to LipidMaps Structure Drawing Tools
`export PATH="\$PATH:./bin/*_lipidmapstools/bin"`;

# Loading Fatty Acids from File
my @FAlist; # A list containing Fatty Acids
open FA_FILE, "./$InputFolderName/FA" or die "Can\'t open \.\/$InputFolderName\/FA \: $!";
while (<FA_FILE>){
	chomp($_);
	@_ = split("\t",$_);
	push (@FAlist, [@_]); 
}
close FA_FILE;

# Sub-routine's

sub uniq {

	my @InArray = @_;
	my @unique = ();
	my %seen   = ();
	foreach my $elem ( @InArray ) {
		next if $seen{ "$elem" }++;
		push @unique, $elem;
	}

	return (@unique);
	# return keys %{{ map { $_ => 1 } @_ }};  
}   

sub DelDuplicates {
	open FILEA, "./$InputFolderName/O1.txt";
	my @list3=<FILEA>;
	close FILEA;
	@list3 = uniq(@list3);
	open FILEB, "> ./$InputFolderName/O1.txt";
	print FILEB @list3;
	close FILEB;
}

sub ModifySMILES {

	$_ = $_[0];

	$_ =~ s/\[C@\@H\]/C/g;			# Removing steriochemistry for chiral carbon
	$_ =~ s/\[C\@H\]/C/g;
	$_ =~ s/\[C\@\@]/C/g;
	$_ =~ s/\[C\@]/C/g;

	$_ =~ s/\\//g;					# Forward and backward slashes (meaning - cis - trans) are removed
	$_ =~ s/\///g;

	$_ =~ s/\[O\-\]/O/g;			# Removing charges for Oxygen and Nitrogen atoms
	$_ =~ s/\[N\+\]/N/g;	

	# Carbon, Nitrogen, Phosphorus, Sulphur, Fluorine and Iodine retain their characters
#	$_ =~ s/O|\[O\-\]/A/g;			# Since character 'O' is not a defined animo acid, Oxyzen is given 'E' character
#	$_ =~ s/Na|\[Na\]|\[Na\+\]/D/g;	# Two character atoms (Na, Br, Cl) are given single character
#	$_ =~ s/Br|\[Br\]/E/g;			
#	$_ =~ s/Cl|\[Cl\]|\[Cl\-\]/G/g;
#	$_ =~ tr/=#/HK/;				# Double and triple bonds are given G and H characters respectively
#	$_ =~ tr/()/LM/;				# Round brackets (meaning - branch opening and closing) are given K and L characters

#	$_ =~ tr/cn/CN/;				# Carbon, Nitrogen in ring structures (represented in small letters) and converted to capital letters
	
#	$_ =~ tr/H//d;					# Hydrogen atoms are removed

#	$_ =~ tr/[+,\-,*,.]//d;			# Charges (+,-), wild character (*) and join character (.) are removed
#	$_ =~ tr/1-9//d;				# Numbers 1-9 are removed
#	$_ =~ tr/[]//d;				# Square brackets are removed; These are used for differentiating special atoms
#	$_ =~ tr/@//d;					# @ (specifies position in with respect to chiral centre) is removed

	return ($_);
}

sub GoToCer { 

	if ($h2==0){ $h_position = "";}
	elsif ($h2==1){ $h_position = "(2OH)";}
	elsif ($h2==2){ $h_position = "(2OH,4OH)";}

	if ($h1==2){ print O1 "$LipidClass\(d$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}
	elsif ($h1==3){ print O1 "$LipidClass\(t$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2$h_position\)\n";}

}

sub Search_UserSMILES {
	my $key;
	foreach $key (keys (%UserSMILES)) { 
		if ($key =~ m/\Q$LipidName\E/){
			# Add2Lists
			$LipidNameCopy =~ s/(\s|;|:|,)/_/g;
			push (@CentroidList, "$LipidNameCopy\tUserDefinedSMILES\t0"); 
			push (@LipidNamesList, "$UserSMILES{$LipidName} $LipidNameCopy");
			if ($LipidName =~ m/_/) {
				@_ = split("_", $LipidNameCopy); push (@ClassList, $_[0]);
			}
			last FOO;
		}
	}
	print "\nCould not find a way to generate structure for this molecule\n\n";
	push @UnsupportedSpecies, $LipidName;
	last FOO;
}

sub Add2Lists {

	# To CentroidList
	$LipidNameCopy =~ s/(\s|;|:|,)/_/g;
	$_ = scalar(@VarList2);
	push (@CentroidList, "$LipidNameCopy\t$VarList2[$centroid]\t$_"); 
	
	# To SMILESlist
	push (@LipidNamesList, "$VarList1[$centroid] $LipidNameCopy"); 

	# To ClassList
	@_ = split("_", $LipidNameCopy);
	push (@ClassList, $_[0]); 

}

sub WriteLists2Files{

	# Write CentroidList to File
	$CentroidListFileName = "CentroidList\.txt";
	open O3, "> ./$InputFolderName/$CentroidListFileName" or die "Unable to open ./$InputFolderName/$CentroidListFileName file";
	print O3 join ("\n", @CentroidList); 
	close O3;
	
	# Write SMILESlist to File
	$SMILESlistFileName = "SMILESlist\.smi";
	open O4, "> ./$InputFolderName/$SMILESlistFileName";
	print O4 join ("\n", @LipidNamesList); 
	close O4; 

	# Write ClassList to File
	@ClassList = uniq(@ClassList); # removing duplicate classes
	open O5, "> ./$InputFolderName/ClassList";
	print O5 join ("\n", @ClassList);
	print O5 "\n";
	close O5;
}

sub CheckForLyso {
	# Finding out of a lipid species is a lyso molecule or not.
	$LysoFlag=0;
	if ($LipidClass =~ m/^TAG/) {$LipidClass =~ s/TAG/TG/;}
	elsif ($LipidClass =~ m/^DAG/) {$LipidClass =~ s/DAG/DG/;}
	elsif ($LipidClass =~ m/^IPC/) {$LipidClass =~ s/IPC/PICer/;}
	elsif ($LipidClass =~ m/^CerPE/) {$LipidClass =~ s/CerPE/PECer/;}
	elsif ($LipidClass =~ m/^MIPC/) {$LipidClass =~ s/^MIPC/MIPCer/;}
	elsif ($LipidClass =~ m/^M\(IP\)2C/) {$LipidClass =~ s/^M\(IP\)2C/MIP2Cer/;}
	elsif ($LipidClass =~ m/^LCB$/) {$LipidClass =~ s/LCB/Cer/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LCBP$/) {$LipidClass =~ s/LCBP/CerP/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPC/) {$LipidClass =~ s/LPC/PC/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPA/) {$LipidClass =~ s/LPA/PA/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPS/) {$LipidClass =~ s/LPS/PS/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPI/) {$LipidClass =~ s/LPI/PI/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LPE/) {$LipidClass =~ s/LPE/PE/; $LysoFlag=1;} # Its a lyso molecule
	elsif ($LipidClass =~ m/^LIPC/) {$LipidClass =~ s/LIPC/PICer/; $LysoFlag=1;} # Its a lyso molecule
}

foreach (@LipidsList){
FOO: {

	$LipidName=$_;
	chomp($LipidName);
	$LipidNameCopy = $LipidName; # making a copy of lipid name to use it later (while generating CentroidList and SMILESlist files )
	
	print "Generating possible structures for $LipidName\n"; # <stdin>; 
	
	# Opening a file to write all isomer possibilities
	open O1, "> ./$InputFolderName/O1.txt"; 
	
	# Breaking lipid name to extract fatty acid chain information
	if ($LipidName =~ m/\s/){
		@_ = split /\s/ , "$LipidName";
		$LipidClass = $_[0]; # print "lipid class = $LipidClass";<stdin>;
		@Chains = split /\// , $_[1];
	}
	else {
		$LipidClass = $LipidName;
	}

	&CheckForLyso;
	
	# Identifying number of c'atoms and d'bonds for (L)P(C|A|S|I|E|G) and DaG and Cer(P) CerPE HexCer
	if ($LipidClass =~ m/\b(PC|PA|PS|PI|PE|PG|DG|Cer|CerP|PICer|HexCer|PECer|MIPCer|MIP2Cer)\b/){

		# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[0]);
		($c1,$d1) = ($_[0], $_[1]);
		if ($LipidClass =~ m/Cer/){ $h1 = $_[2];}
		
		if ($LysoFlag == 1) {
			$c2 = 0;
			$d2 = 0;
			if ($LipidClass =~ m/Cer/){ $h2 = 0;}
		}
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		elsif ($LysoFlag == 0) {
			@_ = split("\:|\;", $Chains[1]);
			($c2,$d2) = ($_[0], $_[1]);
			if ($LipidClass =~ m/Cer/){ $h2 = $_[2];}
		}
		
		foreach (@FAlist){
			if (($_->[0] == $c1) and ($_->[1] == $d1)) {
				if($d1>0){$db_positions1= $_->[2];}else{$db_positions1="";}
				foreach (@FAlist){
					if (($_->[0] == $c2) and ($_->[1] == $d2)) {
						if($d2>0){$db_positions2= $_->[2];}else{$db_positions2="";}
						if ($LipidClass =~ m/DG/){
							print O1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\/0\:0\)\n";
							print O1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\/0\:0\)\n"; # Xchanging SN1 and SN2 positions
						}elsif ($LipidClass =~ m/Cer/){
							&GoToCer($c1,$d1,$h1,$c2,$d2,$h2);
						}else{
							print O1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
							# If its not a Lyso molecule, then exchange SN1 and SN2 positions
							if ($LysoFlag != 1){ 
								# Xchanging SN1 and SN2 positions
								print O1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\)\n"; 
							}
						}
					}
				}
			}
		}
	}
	
	# For TaG's only
	elsif ($LipidClass =~ m/TG/){
		
		# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[0]);	
		($c1,$d1) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[1]);
		($c2,$d2) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[2]);
		($c3,$d3) = ($_[0], $_[1]);

		foreach (@FAlist){
			if (($_->[0] == $c1) and ($_->[1] == $d1)) {
				if($d1>0){$db_positions1= $_->[2];}else{$db_positions1="";}
				foreach (@FAlist){
					if (($_->[0] == $c2) and ($_->[1] == $d2)) {
						if($d2>0){$db_positions2= $_->[2];}else{$db_positions2="";}
						foreach (@FAlist){
							if (($_->[0] == $c3) and ($_->[1] == $d3)) {
							if($d3>0){$db_positions3= $_->[2];}else{$db_positions3="";}
								print O1 "$LipidClass\($c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\)\n";
								print O1 "$LipidClass\($c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\)\n";
								print O1 "$LipidClass\($c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\)\n";
								print O1 "$LipidClass\($c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\)\n";
								print O1 "$LipidClass\($c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\)\n";
								print O1 "$LipidClass\($c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\)\n";
							}
						}
					}
				}
			}
		}
	}									


	# For CL's only
	elsif ($LipidClass =~ m/CL/){

		# Splitting string "chain_1" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[0]);
		($c1,$d1) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[1]);
		($c2,$d2) = ($_[0], $_[1]);
		
		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[2]);
		($c3,$d3) = ($_[0], $_[1]);

		# Splitting string "chain_2" to obtain 1) No. of carbon atoms 2) No. of double bonds
		@_ = split("\:|\;", $Chains[3]);
		($c4,$d4) = ($_[0], $_[1]);

		foreach (@FAlist){
			if (($_->[0] == $c1) and ($_->[1] == $d1)) {
				if($d1>0){$db_positions1= $_->[2];}else{$db_positions1="";}
				foreach (@FAlist){
					if (($_->[0] == $c2) and ($_->[1] == $d2)) {
					if($d2>0){$db_positions2= $_->[2];}else{$db_positions2="";}
					foreach (@FAlist){
							if (($_->[0] == $c3) and ($_->[1] == $d3)) {
							if($d3>0){$db_positions3= $_->[2];}else{$db_positions3="";}
								foreach (@FAlist){
									if (($_->[0] == $c4) and ($_->[1] == $d4)) {
									if($d4>0){$db_positions4= $_->[2];}else{$db_positions4="";}

										print O1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\)\n";

										print O1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\)\n";

										print O1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c4\:$d4$db_positions4\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c4\:$d4$db_positions4\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c3\:$d3$db_positions3\/$c4\:$d4$db_positions4\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\)\n";

										print O1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c3\:$d3$db_positions3\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c1\:$d1$db_positions1\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c2\:$d2$db_positions2\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c3\:$d3$db_positions3\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c2\:$d2$db_positions2\]\, 3\'\-\[$c3\:$d3$db_positions3\/$c1\:$d1$db_positions1\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c1\:$d1$db_positions1\/$c2\:$d2$db_positions2\]\)\n";
										print O1 "$LipidClass\(1\'\-\[$c4\:$d4$db_positions4\/$c3\:$d3$db_positions3\]\, 3\'\-\[$c2\:$d2$db_positions2\/$c1\:$d1$db_positions1\]\)\n";
									}
								}
							}
						}
					}
				}
			}
		}
		
	}
	
	else {
		&Search_UserSMILES;
	}	

	close O1;   
	
	# Program design generates duplicates for Ceramides, PC/PE TAG etc. (TODO : Modify design so that duplicates are not generated)
	# Deleting duplicates...
	&DelDuplicates;
	
	# Generating MOL-MOL file for all possible structures using LipidMaps Structure Drawing Tools
	if ($LipidClass =~ m/(TG|DG)/){ $LipidMapsProgramToUse = "GLStrGen.pl";}
	elsif ($LipidClass =~ m/(Cer|PICer|PECer|HexCer|MIPCer|MIP2Cer)/){ $LipidMapsProgramToUse = "SPStrGen.pl";}
	elsif ($LipidClass =~ m/(PC|PA|PS|PI|PE|PG)/){ $LipidMapsProgramToUse = "GPStrGen.pl";}
	elsif ($LipidClass =~ m/CL/){ $LipidMapsProgramToUse = "CLStrGen.pl";}
	
	`perl ./bin/*_lipidmapstools/bin/$LipidMapsProgramToUse -m AbbrevFileName ./$InputFolderName/O1.txt -r ./$InputFolderName/O1`;
	# `mv O1.sdf ./$InputFolderName/O1.sdf`;
	
	# Converting MOL-MOL structures to SMILES using "OpenBabel : Molecule Conversion package"
	`babel -isdf ./$InputFolderName/O1.sdf -osmi ./$InputFolderName/O1.smi -d -b -r`; # print "check 01.smi file";<stdin>;

	# Finding Centroid structure
		# Changing molecule names (because of limitation in ASCII characters that can be used as Molecule names

		open O1, "./$InputFolderName/O1.smi";
		open O2, "> ./$InputFolderName/O2.smi";
		$i=0;
		@VarList1=(); # Array to store SMILES strings
		@VarList2=(); # Array to store Molecule names
		while (<O1>){
			($VarString1, $VarString2) = split("\t", $_);
			$VarString1 = &ModifySMILES($VarString1);
			push (@VarList1, "$VarString1");
			chomp($VarString2);
			push (@VarList2, "$VarString2");
			print O2 "$VarString1 $i\n";
			$i++;
		}
		close O2;  # print "check 02.smi file";<stdin>;
		close O1;	

	# If there are no possible sructures..
	
	if ($i == 0){ 
		print STDERR "\nUnable to generate structures from Fatty Acid possibilities defined in ./bin/$InputFolderName file\n"; 
		`rm ./$InputFolderName/O1.sdf ./$InputFolderName/O1.smi ./$InputFolderName/O1.txt ./$InputFolderName/O2.smi`;
		&Search_UserSMILES;
	}
	
	# If there is only one possible structure for a given lipid species, take that structure as default centroid structure
	elsif ($i == 1) {$centroid = 0;} 
	
	# If there is only TWO possible structure for a given lipid species, choose fisr one as centroid structure
	elsif ($i == 2) {$centroid = 0;} 
	# If there is only TWO possible structure for a given lipid species, randomly choose one of them as centroid structure
	# elsif ($i == 2) {$centroid = int(rand(2));} 

	# If there are more than TWO possible structures, find centroid  structure using any of the six similarity scoring methods
	else {
	
		if    ($MethodChoice==1){`perl		./bin/*m1.pl	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==2){`python	./bin/*m2.py	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		#elsif ($MethodChoice==3){`perl	./bin/*m3.pl	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==4){`python	./bin/*m4.py	./$InputFolderName/O2.smi					./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==5){`python	./bin/*m5.py	./$InputFolderName/O2.smi	./bin/mat4		./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;}
		elsif ($MethodChoice==6){`perl		./bin/*m6.pl	./$InputFolderName/O2.smi	./$InputFolderName/O2.aln	./$InputFolderName/O2.csv	./$InputFolderName/O2.log`;
						`rm ./$InputFolderName/O2.aln ./$InputFolderName/O2.aln_joined ./$InputFolderName/O2.aln_joined_rearranged`;}

		# print "check csv file";<stdin>;
		# input *.csv file into centroid.pl program
		# get the center molecule and it's SMILES
		# $centroid =`echo "./$InputFolderName/O2.csv" | perl ./bin/centroid2.pl;`; print "$centroid";<stdin>;
		$centroid =`perl ./bin/centroid2.pl ./$InputFolderName/O2.csv`; # print "$centroid";<stdin>;

		# Deleting temporary files..
		`rm ./$InputFolderName/O2.csv ./$InputFolderName/O2.log`;
	
	}

		# Deleting temporary files..
		`rm ./$InputFolderName/O1.sdf ./$InputFolderName/O1.smi ./$InputFolderName/O1.txt ./$InputFolderName/O2.smi`;


	&Add2Lists;
}
}

&WriteLists2Files;

# Calculating pair wise distances between lipid species...

my $CSVfilename = "$InputFolderName"."\.csv"; # print "./$InputFolderName/$CSVfilename";<stdin>;
my $LOGfilename = "$InputFolderName"."\.log";
my $ALNfilename = "$InputFolderName"."\.aln";
	
if    ($MethodChoice==1){`perl		./bin/*m1.pl	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==2){`python	./bin/*m2.py	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
#elsif ($MethodChoice==3){`perl	./bin/*m3.pl	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==4){`python	./bin/*m4.py	./$InputFolderName/$SMILESlistFileName				./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==5){`python	./bin/*m5.py	./$InputFolderName/$SMILESlistFileName	./bin/mat4	./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;}
elsif ($MethodChoice==6){`perl		./bin/*m6.pl	./$InputFolderName/$SMILESlistFileName	./$InputFolderName/O2.aln	./$InputFolderName/$CSVfilename	./$InputFolderName/$LOGfilename`;
					`rm ./$InputFolderName/O2.aln ./$InputFolderName/O2.aln_joined`;`mv ./$InputFolderName/O2.aln_joined_rearranged ./$InputFolderName/$ALNfilename`;}

# Performing Principal Component Analysis 
my $PCAfilename = "$InputFolderName"."\.pca";
my $VariancePlotFileName = "$InputFolderName"."_VariancePlot"."\.pdf";
my @PCAcoordinates = `Rscript ./bin/121012_PCAfromCSV.R ./$InputFolderName/$CSVfilename ./$InputFolderName/$VariancePlotFileName ./$InputFolderName/$PCAfilename`;

chomp($PCAcoordinates[0]); #print "$PCAcoordinates[0] \n";
chomp($PCAcoordinates[1]); #print "$PCAcoordinates[1] \n";
chomp($PCAcoordinates[2]); #print "$PCAcoordinates[2] \n";
chomp($PCAcoordinates[3]); #print "$PCAcoordinates[3] \n $HighestConc \n"; <stdin>;

# Plotting programs
my @MultiPlotOutput=();

foreach $TempFileName (@FileNameList){
	@MultiPlotOutput = `Rscript ./bin/130115_PCA_plot.R "./$InputFolderName/$TempFileName" "./$InputFolderName/$InputFolderName\.pca" "./$InputFolderName/ClassList" "./$InputFolderName/CentroidList.txt" "$PCAcoordinates[0]" "$PCAcoordinates[1]" "$PCAcoordinates[2]" "$PCAcoordinates[3]" "$HighestConc" "./$InputFolderName/$TempFileName"`;
	# print @MultiPlotOutput;
}

# Combined Lipidome Plot
my $SVGfilename_2D = "$InputFolderName"."_2D\.svg";
my $SVGfilename_3D = "$InputFolderName"."_3D\.svg";
`Rscript ./bin/130124_CombinedLipidomePlot.R "./$InputFolderName/$SVGfilename_2D" "./$InputFolderName/$PCAfilename" "./$InputFolderName/ClassList" "./$InputFolderName/$SVGfilename_3D"`;

# Writing un-analyzed lipid species to <stdout>
print "\nUnused lipid abbreviations\n\n";
print join ("\n", @UnsupportedSpecies);
print "\n\n";

# Calculating Hausdroff distance and Average distance from PCA Coordinates
my $FileNameListRef = \@FileNameList;
$FileNameList = join("XXXX", @FileNameList);
`perl ./bin/130123_MaxDist_and_AvgDist_from_PC1_PC2_PC3.pl $PCAfilename $FileNameList $InputFolderName`;

# Calculating Hausdroff distance and Average distance from Similarity Scores
$FileNameList = join("XXXX", @FileNameList);
`perl ./bin/130117_MaxDist_and_AvgDist_from_SimilarityScore.pl $CSVfilename $FileNameList $InputFolderName`;

# Creating Dendrograms for Pairwise 1) Overlap Score (from Pairwise Distances - Similarity scoring methods) 2) Hausdroff distance (from PCA Coordinates) 3) Average distance
`Rscript ./bin/130117_Hdist2Dendro.R $InputFolderName "./$InputFolderName/WithConcMaxdist_from_PC1_PC2_PC3.csv" "./$InputFolderName/WithOutConcMaxdist_from_PC1_PC2_PC3.csv" "./$InputFolderName/WithConcAvgDist_from_PC1_PC2_PC3.csv" "./$InputFolderName/WithOutAvgDist_from_PC1_PC2_PC3.csv" "./$InputFolderName/WithConcMaxdist_from_SimilarityScore.csv" "./$InputFolderName/WithOutConcMaxdist_from_SimilarityScore.csv" "./$InputFolderName/WithConcAvgDist_from_SimilarityScore.csv" "./$InputFolderName/WithOutAvgDist_from_SimilarityScore.csv"`;


__END__
