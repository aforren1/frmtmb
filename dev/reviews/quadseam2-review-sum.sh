#!/usr/bin/env bash
# Review of lane quadseam2: the summary lines the review quotes, generated
# from the logs in dev/reviews/quadseam2-review-log/.
cd /c/Users/adf44/source/r/frmtmb-wt-quadseam2/dev/reviews/quadseam2-review-log
echo "== ident (attack 1, 2)"
for k in FIT SE; do
  diff <(grep "^$k " ident-base.txt) <(grep "^$k " ident-lib.txt) > /dev/null &&
    echo "$k lines identical: $(grep -c "^$k " ident-lib.txt)" || echo "$k DIFFER"
done
echo "== esc2 (attack 4)"
join <(grep ^ROW esc2-base.txt | awk '{print $2"/"$3, $5, $7, $9, $11}' | sort) \
     <(grep ^ROW esc2-lib.txt | awk '{print $2"/"$3, $5, $7, $9, $11}' | sort) |
  awk '{d=$4-$8; if (d<0) d=-d; n++; if (d>1e-8) {k++; printf "%-12s base rule %s ll %s alpha %s | lane rule %s ll %s alpha %s | lost %.4f\n",$1,$3,$4,$5,$7,$8,$9,$4-$8}} END{print "fits",n,"lane below base by >1e-8:",k}'
echo "== limit (attack 3)"
grep ^ROW limit-lib.txt | awk -F'\t' '{k=$2"/"$3; v=$7; sub(/ maxse.*/,"",v); gsub(/frmtmb_warning,?/,"",v); gsub(/boundary\+separation/,"separation+boundary",v); if($4=="LAP") L[k]=v; else Q[k]=v} END{n=0; d=0; for(k in L){n++; if(L[k]!=Q[k]) {d++; print "DIFF",k}} print "lane quad vs Laplace verdicts: pairs",n,"differ",d}'
join -t $'\t' <(grep ^ROW limit-base.txt | awk -F'\t' '$4=="QUAD"{print $2"/"$3"\t"$7}' | sort) \
  <(grep ^ROW limit-lib.txt | awk -F'\t' '$4=="QUAD"{print $2"/"$3"\t"$7}' | sort) |
  awk -F'\t' '{a=$2; b=$3; sub(/ maxse.*/,"",a); sub(/ maxse.*/,"",b); if (a!=b) {c++; if (a ~ /frmtmb_warning/ && b !~ /frmtmb_warning/) g++}} END{print "quad verdicts base vs lane differ:",c,"of which base had the extra frmtmb_warning:",g}'
echo "== touch (attack 2)"
for s in "" _b; do
  join -t $'\t' <(grep ^ROW touch-base$s.txt | awk -F'\t' '{print $2"/"$3"/"$4"\t"$5"\t"$9"\t"$10"\t"$11"\t"$12"\t"$8}' | sort) \
    <(grep ^ROW touch-lib$s.txt | awk -F'\t' '{print $2"/"$3"/"$4"\t"$5"\t"$9"\t"$10"\t"$11"\t"$12"\t"$8}' | sort) |
    awk -F'\t' -v s="$s" '{n++; if($2==$8) same++; if($5!=$11||$6!=$12) v++; if($13!="NA") e++; bf+=$3; bv+=$4; lf+=$9; lv+=$10} END{printf "run%s fits %d logLik bitwise same %d verdicts differ %d restarts kept %d | secs base fit %.1f vcov %.1f lane fit %.1f vcov %.1f\n",s,n,same,v,e,bf,bv,lf,lv}'
done
echo "== passes (attack 1)"
join <(grep ^ROW passes-base.txt | awk '{print $2"/"$3, $5, $11, $15}' | sort) \
     <(grep ^ROW passes-lib.txt | awk '{print $2"/"$3, $5, $11, $13, $15}' | sort) |
  awk '{printf "%-11s base starts %s ll %s fit %ss | lane starts %s ll %s free %s fit %ss | lane-free %.3g\n",$1,$2,$3,$4,$5,$6,$7,$8,$6-$7}'
echo "== latres (attack 5)"
tok() {
  grep -v "^lib\|^END" "$1" | awk '{s=""; for (i=1;i<NF;i++) if ($i=="latres") s=s" "$(i+1); print s}'
}
paste -d"|" <(tok latres-base.txt) <(tok latres-lib.txt) \
  <(grep -v "^lib\|^END" latres-lib.txt | cut -c1-60) |
  awk -F"|" '$1!=$2{print "CHANGED base", $1, "lane", $2, ":", $3} END{print "latres rows compared:", NR}'
grep "SEED\|gate" fewbern-lib.txt
