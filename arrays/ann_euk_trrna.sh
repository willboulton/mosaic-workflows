swrap --partition compute-64-512,compute-24-96 --mem 32G -n 1 -c 16 --wrap "conda activate trnascan_se; cd /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/annotations/euk; ls repeatmasker | parallel -j 16 \"barrnap --kingdom euk --quiet repeatmasker/{}/{}.fna.masked | tail -n +2 > rrna/{}/{}.gff\""

swrap --partition compute-64-512,compute-24-96 --mem 64G -n 1 -c 16 --wrap "conda activate trnascan_se;cd /gpfs/data/mock_lab/will/mosaic_tmp/mag_catalogue/annotations/euk;ls repeatmasker | parallel -j 16 \"tRNAscan-SE -E --thread 1 --quiet --gff trna/{}/{}.gff repeatmasker/{}/{}.fna.masked\""
