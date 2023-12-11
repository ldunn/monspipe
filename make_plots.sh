#!/usr/bin/env bash

psrj=$1
date=$2

ln -s ${psrj}_${date}.ar ${date}.ar
ln -s ${psrj}_${date}.clean ${date}.clean

pav -GTd ${date}.ar -g ${psrj}_${date}_GT_raw.png/PNG
pav -GTd ${date}.clean -g ${psrj}_${date}_GT_clean.png/PNG

pav -FYd ${date}.ar -g ${psrj}_${date}_FY_raw.png/PNG
pav -FYd ${date}.clean -g ${psrj}_${date}_FY_clean.png/PNG

pav -DFTp ${date}.ar -g ${psrj}_${date}_DFTp_raw.png/PNG
pav -DFTp ${date}.clean -g  ${psrj}_${date}_DFTp_clean.png/PNG

rm ${date}.ar
rm ${date}.clean
