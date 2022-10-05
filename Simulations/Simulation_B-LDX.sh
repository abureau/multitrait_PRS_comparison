source s-ldxr-0.3-beta/s-ldxr_env/bin/activate
cd s-ldxr-0.3-beta
python s-ldxr.py \
    --gcor .../S-LDXR/ss_SKZ.txt.gz \
           .../S-LDXR/ss_BIP.txt.gz \
    --ref-ld-chr .../baseline-LD-X/ldscore/EAS_EUR_allelic_chr \
    --w-ld-chr .../baseline-LD-X/ldscore/EAS_EUR_allelic_chr \
    --frqfile .../maf/EUR/1000G.EUR.QC. \
              .../maf/EUR/1000G.EUR.QC. \
    --annot .../baseline-LD-X/annotations/baseline-LD-X. \
    --apply-shrinkage 0.5 \
    --save-pseudo-coef \
    --out .../S-LDXR/out.txt
deactivate