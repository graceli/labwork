cd kk/
new
cd labtalk1_11232006/
new
open 300_50ns_pmf.bmp
open grace_2006.ppt 
open -a keynote grace_2006.ppt 
new
cd 
new
cd Desktop/
new
open .
new
s
ssh t2.research.sickkids.ca
exit
lat
ssh lattice.westgrid.ca
ssh t2.research.sickkids.ca
c
cd labwork/
cd ../
cd github/papers/
new
cd abeta_full/
enw
new
mkdir ~/bin/textmate
cp ~/Downloads/Completion/BibDeskTMCompletions ~/bin/textmate
new
mate ab-body.tex
mate ab-body.tex
open -a ab.bib 
open -a bibdesk ab.bib 
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
mate ab-body.tex
~/bin/textmate/BibDeskTMCompletions 
file ~/bin/textmate/BibDeskTMCompletions 
~/bin/textmate/BibDeskTMCompletions | Less
~/bin/textmate/BibDeskTMCompletions | less
ls
cd ../
~/bin/textmate/BibDeskTMCompletions | less
cd
~/bin/textmate/BibDeskTMCompletions | less
env  | grep TM
echo $TM_SUPPORT_PATH
cd /
find . -name exit_codes.rb
~/bin/textmate/BibDeskTMCompletions -h
vi /System/Library/Frameworks/Ruby.framework/Versions/1.8/usr/lib/ruby/1.8/pathname.rb
ssh tera.research.sickkids.ca
cd labwork/
ls
st
cat
s
ssh t2.research.sickkids.ca
c
alias
alias
c
vi ~/.bashrc 
c
re
new
find . -name *biophysical*
open ./KLV_paper/starting_pieces/biophysical_poster_printed.pdf
find . -name *seminar*
man zip
zip
unzip
unzip -c archives/biophysical_2011.zip 
unzip -l archives/biophysical_2011.zip 
man unzip
unzip -l archives/biophysical_2011.zip 
unzip -p archives/biophysical_2011.zip biophysical_2011/biophysical_2011.key
unzip -x archives/biophysical_2011.zip biophysical_2011/biophysical_2011.key
new
open biophysical_2011/biophysical_2011.key 
find . -name *seminar*
unzip -l ./archives/seminar2012.zip
unzip -x ./archives/seminar2012.zip seminar2012/student_seminar_2012_final.key
open seminar2012/student_seminar_2012_final.key
find . -name *hpcs*
cd ./work/2010/may/May_2010/hpcs/hpcs_poster
ls
cd hpcs_new_renders/
new
open *.tga
new
open *.png
cd labwork/
s
s
cd labwork/
new
find . -name *inositol*
find . -name *inositol*.sh
cat ./abeta_analysis/scripts/inositol_binding_extract.sh
cat ./klv_analysis/setup/add_inositol.sh 
find . -name *ions.sh
find . -name *ions*.sh
find . -name *ions*
re
find . -name *add*
find . -name *inositol*
find . -name *inositol* | grep sh
cd ~/labwork/
ls
cd abeta_analysis/
ls
ls *.sh
cd scripts/
ls
cd setup/
ls
cat setup.sh 
mate setup.sh 
pwd
re
new
cd ab_paper/
new
cd resources/
ls
cd ../
new
mkdir runs
cd runs/
mkdir epi
cd epi/
scp -rp colosse:/rap/uix-840-ac/grace/abeta/42/glucose/sys1/em.gro .
scp -rp colosse.clumeq.ca:/rap/uix-840-ac/grace/abeta/42/glucose/sys1/em.gro .
vmd em.gro 
c
c
s
c
c
lat
alias
ssh lattice
ssh lat
exit
re
cd ab_paper/
new
cd volmap
new
ls *.png
open *.png
open *.png
cd ../
ls
new
find . -name *.png
find . -name "*.png"
open ./analysis_march_2012/volmap_xtcs/scyllo_64_test_hardplastic.dat.png 
cd ./analysis_march_2012/volmap_xtcs/
new
open *.png
cd ../
cd ../
new
ls check_tpr_is_crystl/
rm -rf check_tpr_is_crystl/
new
open *.png
ls *.png
rm chiro.png scyllo.png glycerol.png 
new
alias 
alias | grep rm
new
cat plot.py 
new
ptdump  test.h5 
ptdump -v test.h5 
ptdump -v analysis.h5 
new
vi ~/.bashrc 
. ~/.bash_aliases 
rm test.h5 
new
mkdir figures
mv nmr_structure.png figures/
new
ls glycerol/
new
mkdir resources
mv *.key resources/
mv committee_meeting_oct22_2010 resources/
new
mv *.ppt resources/
new
open abeta_literature_notes.pages 
man ruby
new
rm abeta_literature_notes.pages 
rm Backup\ of\ abeta_literature_notes.pages 
new
new
cd runs/
new
cd epi/
new
cat em.gro | grep -v INS
cat em.gro | grep -v INS | grep -v CL | grep -v NA
cat em.gro | grep -v INS | grep -v CL | grep -v NA | grep -v SOL 
cat em.gro | grep -v GLCA | grep -v CL | grep -v NA | grep -v SOL 
cat em.gro | grep -v GLCA | grep -v CL | grep -v NA | grep -v SOL  > protein.gro
vi protein.gro 
cat protein.gro | wc -l
vi protein.gro 
new
new
vmd protein.gro 
new
scp -rp colosse.clumeq.ca:/rap/uix-840-ac/grace/abeta/42/inositol/epi_setup/15/epi/1/ions_counter.gro .
vmd ions_counter.gro 
scp -rp colosse.clumeq.ca:/rap/uix-840-ac/grace/abeta/42/inositol/epi_setup/64/epi/1/ions_counter.gro 64.gro
vmd 64.gro 
scp -rp colosse.clumeq.ca:/rap/uix-840-ac/grace/abeta/42/inositol/epi_setup/64/epi/2/ions_counter.gro 64.gro
vmd ions_counter.gro 
re
cd ab_paper/
new
cd runs/epi/
new
vi /Users/grace/labwork/abeta_analysis/scripts/setup/setup.sh
vi /Users/grace/labwork/abeta_analysis/scripts/setup/setup.sh
/Users/grace/labwork/abeta_analysis/scripts/setup/setup.sh protein.gro 
new
ls 64/
ls 64/epi/
vi /Users/grace/labwork/abeta_analysis/scripts/setup/setup.sh
new
rm -rf 64 15/
/Users/grace/labwork/abeta_analysis/scripts/setup/setup.sh protein.gro 
new
scp -rp colosse:/rap/uix-840-ac/grace/abeta/42/glucose/abeta42_glucose.top .
scp -rp colosse.clumeq.ca:/rap/uix-840-ac/grace/abeta/42/glucose/abeta42_glucose.top .
cat abeta42_glucose.top 
vi abeta42_glucose.top 
vi abeta42_glucose.top 
new
rm \#*
new
rm -rf 15
new
cd ../
new
scp -rp epi colosse.clumeq.ca:/rap/uix-840-ac/grace/abeta/42/inositol/epi_setup
history 
new
cd epi/
new
c
c
c
c
c
c
c
s
cd Desktop/
open .
/Users/grace/Desktop/me.bib 
c
cd github/
ls
new
cd papers/
ls
new
mkdir inos2
new
cd inos2/
new
cp ../abeta_full/* .
cp -rp ../abeta_full/* .
new
mate *.tex
new
grep * "ab'
grep * "ab"
grep "ab" *
grep "ab" *.sh
grep -n "ab" *.sh
ls *.sh
mate makepdf.sh 
new
mate makepdf.sh 
cp ../abeta_full/makepdf.sh .
mate makepdf.sh 
./makepdf.sh ab
st
cd ../
git add .
git commit -am "Added the first latex skeleton for inos2"
git push origin master
new
cd inos2/
new
open klvffae.bib 
open -a bibdesk klvffae.bib 
mv ~/Desktop/klvffae.bib .
open -a bibdesk klvffae.bib 
./makepdf.sh ab
new
cat makepdf.sh 
rm ab.bib 
new
./makepdf.sh ab
new
mv klvffae.bib ab.bib
vi makepdf.sh 
./makepdf.sh ab
new
rm -rf build/
./makepdf.sh ab
new
cd build/
ls
new
less ab.bib b
less ab.bib
cd ../
new
mate ab.bib 
./makepdf.sh ab
cd github/papers/inos2/
new
rm -rf build/
mate ../inos2/
new
./makepdf.sh 
./makepdf.sh ab
new
cd build/
new
less ab.aux 
pwd
cd ../
./makepdf.sh ab
./makepdf.sh ab
./makepdf.sh ab
cd github/
cd papers/inos2/
ls
new
open -a bibdesk ab.bib 
mate ../inos2/
./makepdf.sh ab
new
cat makepdf.sh 
grep klvffae *.tex
cd ../
cd abeta_full/
ls
new
mate ab-body.tex 
c
alias
vi ~/.bash_aliases 
. ~/.bash_aliases 
dbox
cd Documents/work/in_progress/
cd inositol_paper/
open -a mou scratch_writing.md 
open -a mou
open -a mou.app
open -a Mou scratch_writing.md 
open -a Mou
open scratch_writing.md 
open -a keynote
new
mate scratch_storyboarding.md
mate scratch_discussion.md
new
cd inos
ls *.enl
ls
cd inositol_
cd inositol_bib/
ls
new
open klv_paper_clean3.enl
c
alias
vi ~/.bashrc 
vi ~/.bash_profile
vi ~/.bash_aliases 
vi ~/.bashrc 
. ~/.bash_aliases 
alias
is
vi ~/.bash_aliases 
vi ~/.bash_aliases 
cd labwork/
st
alias
cd ..
cd systems/
new
st
git pull
cd ~/labwork/
st
git commit -am "epi setup"
git push origin master
cd github/papers/
new
st
cd abeta_full/
st
git status
mate ../abeta_full/
new
./makepdf.sh 
./makepdf.sh 
./makepdf.sh 
./makepdf.sh 
./makepdf.sh 
./makepdf.sh 
./makepdf.sh 
./makepdf.sh 
st
git commit -am "Added a paragraph of discussion of the sig. result of binding to KLVFFAE face"
git push origin master
open file://localhost/Users/grace/Documents/Papers2/Articles/2010/Wu/Wu_2010-1.pdf
open /Users/grace/Documents/Papers2/Articles/2010/Wu/Wu_2010-1.pdf
open /Users/grace/Documents/Papers2/Articles/2010/Wu/Wu_2010.pdf 
open /Users/grace/Documents/Papers2/Articles/2010/Wu/Wu_2010=2.pdf 
open /Users/grace/Documents/Papers2/Articles/2010/Wu/Wu_2010-2.pdf 
open file://localhost/Users/grace/Documents/Papers2/Articles/2010/Wu/Wu_2010-1.pdf
open /Users/grace/Documents/Papers2/Articles/2010/Wu/
st
git status
git commit -am "mod ab.bib"
git push origin master
st
git status
mate ../abeta_full/
git commit -am "Added some literature notes"
git push origin master
git commit -am "Added more comments on the hypothesis coming out of this work"
git push origin master
s
s
cd github/
cd papers/
git pull
ls
mate sh3/
git status
git commit -am "Added comments on why sh3-denatured more compact than sh3-nondenatured"
git push origin master
new
mate sh3/
cd sh3/
new
./makepdf.sh sh3
st
git commit -am "Added some notes from Sh3 meeting 1. Not yet complete"
git push origin master
st
cd github/
ls
new
cd papers/
new
cd inos2/
git status
git pull origin master
st
git log
st
git push origin master
st
lat
ssh lattice
ssh lattice.westgrid.ca
exit
c
exit
s
exit
ssh goblin.sharcnet.ca
ls
