for ((i=100;i;i--));
    do ./merge ./Run$i/AvaTot.tsv AvaTot.tsv
    ./merge ./Run$i/AvaPEle.tsv AvaPEle.tsv
    ./merge ./Run$i/EleEner.tsv EleEner.tsv
    ./merge ./Run$i/ElePPhot.tsv ElePPhot.tsv
    ./merge2 ./Run$i/EEnerPos.tsv EEnerPos.tsv
    ./merge ./Run$i/ElePos.tsv ElePos.tsv;
done
