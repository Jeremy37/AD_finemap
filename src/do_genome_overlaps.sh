#!/bin/bash

F=$1
NAME=$2

zcat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3}' > overlaps/$NAME.overlap_input.bed
zcat $F | sed '1d' | awk 'BEGIN{OFS="\t"}{print "chr"$1,$2-1,$2,$3}' > overlaps/$NAME.overlap_input.chr.bed

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.txt --output overlaps/roadmap.dnase -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.brain.txt --output overlaps/roadmap.dnase.brain -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.dnase.blood_and_immune.txt --output overlaps/roadmap.dnase.blood_and_immune -vv

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/fantomEnh.bedfile.txt --output overlaps/FantomEnh -vv

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.txt --output overlaps/RoadmapEnh -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.brain.txt --output overlaps/RoadmapEnh.brain -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/roadmap.bedfilelist.RoadmapEnh.blood_and_immune.txt --output overlaps/RoadmapEnh.blood_and_immune -vv

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/ipsc.atac.bedfile.txt --output overlaps/ipsc.atac -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/npc.atac.bedfile.txt --output overlaps/npc.atac -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/neuron.atac.bedfile.txt --output overlaps/neuron.atac -vv
multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/ineuron.atac.bedfile.txt --output overlaps/ineuron.atac -vv

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.chr.bed --bedfilelist overlaps/ipsMacrophage.atac.bedfile.txt --output overlaps/ipsMacrophage.atac -vv

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.bed --bedfilelist overlaps/ipsSensoryNeuron.atac.bedfile.txt --output overlaps/ipsSensoryNeuron.atac -vv

multi_bedintersect.py --snpbed overlaps/$NAME.overlap_input.bed --bedfilelist overlaps/microglia.atac.bedfile.txt --output overlaps/microglia.atac -vv


# Merge together all the results
paste <(zcat $F | head -n 1 | cut -f 3) \
      <(echo -e "iPSC\tmicroglia\tipsMacrophage\tNPC\tNeuron\tiNeuron\tipsSensNeuron\tDNase\tBrain DNase\tBlood & Immune DNase\tFantom Enh\tRoadmapEnh\tBrain Enh\tBlood & Immune Enh") \
      <(echo -e "DNase overlaps\tBrain DNase overlaps\tBlood & Immune DNase overlaps\tRoadmapEnh overlaps\tBrain Enh overlaps\tBlood & Immune Enh overlaps") \
      > annotated/$NAME.overlaps.txt
paste <(zcat $F | sed '1d' | cut -f 3) \
      <(sed '1d' overlaps/ipsc.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/microglia.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/ipsMacrophage.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/npc.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/neuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/ineuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/ipsSensoryNeuron.atac.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.brain.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.blood_and_immune.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/FantomEnh.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.brain.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/RoadmapEnh.blood_and_immune.overlap.summary.txt | cut -f 4) \
      <(sed '1d' overlaps/roadmap.dnase.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/roadmap.dnase.brain.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/roadmap.dnase.blood_and_immune.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.brain.overlap.summary.txt | cut -f 5) \
      <(sed '1d' overlaps/RoadmapEnh.blood_and_immune.overlap.summary.txt | cut -f 5) \
      >> annotated/$NAME.overlaps.txt
