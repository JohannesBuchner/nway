python3 genrandom_geometric.py || exit 1

topcat -stilts tcopy in=randomcatX.csv ifmt=CSV out=randomcatX.fits || exit 1
nway-write-header.py randomcatX.fits CHANDRA 0.010 || exit 1

topcat -stilts tcopy in=randomcatO.csv ifmt=CSV out=randomcatO.fits || exit 1
nway-write-header.py randomcatO.fits OPT 0.010 || exit 1

nway.py --radius=10.0 randomcatX.fits :pos_err randomcatO.fits 0.1 --out=random_circtest.fits --min-prob=0.01 || exit 1

nway.py --radius=10.0 randomcatX.fits :a:b randomcatO.fits 0.1 --out=random_asymtest.fits --min-prob=0.01 || exit 1

nway.py --radius=10.0 randomcatX.fits :a:b:phi randomcatO.fits 0.1 --out=random_elltest.fits --min-prob=0.01 || exit 1


echo topcat -stilts plot2sky \
   xpix=1081 ypix=548 \
   crowd=0.9998301109057076 \
   clon=150.0338559501912 clat=2.040905749390326 radius=0.003505914272661 \
   auxmin=0 auxmax=0.1 \
   auxvisible=true auxlabel=p_i auxcrowd=0.9998301109057076 \
   legend=true \
   in=random_elltest.fits \
    leglabel='1: All' \
   layer_1=Mark \
      lon_1=CHANDRA_RA lat_1=CHANDRA_DEC \
      shading_1=auto \
   layer_2=SkyEllipse \
      lon_2=CHANDRA_RA lat_2=CHANDRA_DEC ra_2=CHANDRA_a/60/60 \
       rb_2=CHANDRA_b/60/60 posang_2=CHANDRA_phi \
      shading_2=auto \
   layer_3=Mark \
      lon_3=OPT_RA lat_3=OPT_DEC aux_3=p_i \
      shading_3=aux size_3=3 \


nway-explain.py random_circtest.fits 95 || exit 1
nway-explain.py random_asymtest.fits 95 || exit 1
nway-explain.py random_elltest.fits 95 || exit 1
