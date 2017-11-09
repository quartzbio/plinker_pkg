dataset downloaded from http://www.cog-genomics.org/plink/1.9/resources#teach
and converted to bed using:
```
plink --file example/extra --make-bed
```
Then the plink.fam file has been tweaked to include more cases:
 - sex/no sex
 - founder/non-founders
 - no FID
 - duplicated FID
 - parent ID not in fam file


Details:

 - #1 CH18526 NA18526 - set sex to 0 (missing)
 - #2 set Family name (FID) to 0
 - #3 set Family name (FID) to 0
 - #4 set Family name to BOND, and IID to JAMES
 - #5 set Family name to BOND, and IID to ARLETTE
 - #6 set Family name to BOND, IID to VAGA and father ID to JAMES (non-founder)
 - #7 set Family name to BOND, IID to PUDI and mother ID to ARLETTE (non-founder)
 - #8 set Family name to BOND, IID to RE and mother ID to FURI [not in fam] (non-founder??)


 The original file is saved in plink.fam.original.
