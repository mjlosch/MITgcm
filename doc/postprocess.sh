#!/bin/bash
newbib=mr_new.bib
#betterbib manual_references.bib mr.bib
\cp -f mr.bib tmp.bib
# replace umlaute and other special characters
sed -i .bak 's/ö/{\\"o}/g' tmp.bib
sed -i .bak 's/ü/{\\"u}/g' tmp.bib
sed -i .bak 's/ä/{\\"a}/g' tmp.bib
sed -i .bak 's/ç/{\\v{c}}/g' tmp.bib
sed -i .bak "s/á/{\\\'a}/g" tmp.bib
sed -i .bak "s/à/{\\\`a}/g" tmp.bib
sed -i .bak "s/é/{\\\'e}/g" tmp.bib
sed -i .bak "s/è/{\\\`e}/g" tmp.bib

# replace "-" by "--" only in pages fields
# but only if is not preceeded or followed by a "-",
# the following lines do not do it right, therefore first
# replace all "--" by"-"
sed -i .bak '/pages/s/--/-/g' tmp.bib
sed -i .bak '/pages/s/-/--/g' tmp.bib
#
sed -i .bak '/number/s/--/-/g' tmp.bib
sed -i .bak '/number/s/-/--/g' tmp.bib

# replace journals by abbreviations (should be a python script!!!!)
grep "@string" manual_references.bib \
    | sed 's/~/ /g' | sed 's/{{/{/g' | sed 's/}}}/}}/g' > ${newbib}

#    sed -i .bak ":$1:s:.*: journal = $2,:g" tmp.bib
#    sed -n -i .bak "/$1.*journal/s/.*/ journal = $2,/g" tmp.bib
#    sed -i .bak "/$1/s/.*/ journal = $2,/g" tmp.bib
replace_jrnl () {
    sed -i .bak "/^ journal = {$1}/s/.*/ journal = $2,/g" tmp.bib
}
add_and_replace_jrnl () {
    sed -i .bak "/$1/s/.*/ journal = $2,/g" tmp.bib
    echo "@string{$2 = {$1}}" >> ${newbib}
}

listOfAbbrv=`awk -F{ '{print $2}' ${newbib} | awk -F= '{print $1}'`
for abbrv in $listOfAbbrv; do
    jrnl=`awk -F= /"{$abbrv "/'{print $NF}' ${newbib} | sed 's/{//g' | sed 's/}//g' | sed 's/^[ \t]*//' | sed 's/\\ / /g'`
    #    echo "$abbrv = {$jrnl}"
    replace_jrnl "${jrnl}" ${abbrv}
done
# and now the special cases:
replace_jrnl "Ocean Modell." om
replace_jrnl "Q.J Royal Met. Soc." qjrms
replace_jrnl "Phys. Fluids" pf
replace_jrnl "J.~Phys.~Oceanogr." jpo
replace_jrnl "Philosophical Transactions of the Royal Society A: Mathematical, Physical and Engineering Sciences" ptrsl
replace_jrnl "Oceanog." ocn
replace_jrnl "J. Climate" jclim
replace_jrnl "Climate Dyn." climd
replace_jrnl "J. Mar. Syst." jms
replace_jrnl "SIAM J. Sci. Comput." siscicom
replace_jrnl "Tellus A: Dynamic Meteorology and Oceanography" tel

# new abbreviations
add_and_replace_jrnl "The Cryosphere" tc
add_and_replace_jrnl "J. Atmos. Oceanic Technol." jtech
add_and_replace_jrnl "Geosci. Model Dev." gmd
add_and_replace_jrnl "Bull. Amer. Meteor. Soc." bams
add_and_replace_jrnl "Optim. Method. Softw." oms
add_and_replace_jrnl "ACM Trans. Math. Softw." atms
add_and_replace_jrnl "Global Biogeochem. Cycles" gbc
add_and_replace_jrnl "Rev. Geophys." rgp
#add_and_replace_jrnl "Tellus A: Dynamic Meteorology and Oceanography" tela
add_and_replace_jrnl "Deep Sea Res. Part II" dsrtwo
add_and_replace_jrnl "Deep Sea Research Part A. Oceanographic Research Papers" dsra

cat tmp.bib >> ${newbib}
\rm -f tmp.bib*

# comment
# needs manual fixing:
# - betterbib converts hill:99 converts to the content of hoe:99
# - in hill:04 all co-authors are stripped
# still no solution for proper (automatic) capitalization, I may need to do
# this manually
