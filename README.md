# Planetary Magnetic Field Line Tracer

The mag_trace function propogates the magnetic field line from user defined starting position(s), using the planetary internal and external field model functions (also provided by the user). Additional optional inputs can be provided to change default settings: 
    Altitude for field line termination,
    Planet/body radii of triaxial ellipsoid,
    Max propogation distance,
    Max radius of field line from planet/body,
    Min radius for applying external magnetic field model

The vectorized function allows single or multiple input positions. Sets of field line points are output as a cell array, with a field line array for each input position. The current configuration uses the native ode45 integrator, leveraging the variable step size propogation for faster computation, as well as the built-in termination functionality. Other integrators are also being considered for easier translation to other platforms (python, IDL, etc). 

Authors: M. Brennan C. Lawler, and R.J. Wilson

Code Citation DOI: 10.5281/zenodo.10018966 (For all versions.)
Each version release also has its own DOI, click the link above to get to the DOI of specific versions.

This is part of a community code project: Magnetospheres of the Outer Planets Group Community Code and is formulated to work seamlessly with the Jupiter internal and external magnetic field models as outlined in the January 2023 in Space Science Reviews paper:

Wilson, R.J., Vogt, M.F., Provan, G. et al. Internal and External Jovian Magnetic Fields: Community Code to Serve the Magnetospheres of the Outer Planets Community. Space Sci Rev 219, 15 (2023). https://doi.org/10.1007/s11214-023-00961-3

Also see [Magnetospheres of the Outer Planets Group Community Code](https://lasp.colorado.edu/mop/missions/juno/community-code/) for source code links and other references.

# Examples:
    Example 1 Input:
        x = [0;0;0;0;  0];
        y = [5;6;7;8;199];
        z = [0;0;0;0;  0];
        int_field_function = @(x,y,z) (jovian_jrm33_order13_internal_xyz(x,y,z));
        ext_field_function = @(x,y,z) (con2020_model_xyz('analytic',x,y,z));
        [maglinen, maglines] = mag_trace(x,y,z, int_field_function, ext_field_function);
        for n=1:length(maglinen)
            magfpn(n,:) = maglinen{n}(end,:);
            magfps(n,:) = maglines{n}(end,:);
        end
    Example 1 Output:
        magfpn =
        -0.212619628805482         0.346419672517159         0.860448148699362
        -0.227286988204301         0.321669827966839         0.865560434844016
        -0.238198306790462         0.302900806157292         0.868912449725239
        -0.246527424543801         0.288501635104101         0.871164750018913
        -0.431503710396772       -0.0179081056493706         0.849546048892076
        magfps =
        0.0950727966258886         0.340672575791845          -0.8806127755025
        0.0955572881083896         0.306168780506791        -0.891588313656359
        0.0960974912093074          0.28097992615566        -0.898766942642357
        0.0965802754745113         0.262170298626307        -0.903682127540211
        0.0225115478964109       -0.0580771212356769        -0.938914288889963

        maglinen{1} = 
                             0                         5                         0
         -2.04176581735287e-05          4.99999323924756      0.000126191197748637
         -4.08353595166567e-05          4.99998646850128      0.000252381852690593
           -6.125310390317e-05          4.99997968776125      0.000378571964058066
         -8.16708912068563e-05          4.99997289702759      0.000504761531083277
         -0.000183760467065564          4.99993879345838       0.00113570117420451
         -0.000285851096922687          4.99990444006399       0.00176662709361397
         -0.000387942765003547          4.99986983685758       0.00239753919336242
         -0.000490035455534521          4.99983498385241       0.00302843737751371
          -0.00100051369294183          4.99965697231366       0.00618271620748386
          -0.00151101512617662          4.99947271783271       0.00933663276468375
          -0.00202153778478639          4.99928222214286        0.0124901750721737
          -0.00253207969905454          4.99908548702021        0.0156433311616372
          -0.00508500922852381           4.9980082850458        0.0314029001210928
          -0.00763812515488526          4.99677538819522         0.047151023577404
           -0.0101911822946513          4.99538706616396        0.0628862163960519
            -0.012743936162928          4.99384361482233        0.0786070002439256
           -0.0254946787456282          4.98381146286566         0.156943692414671
           -0.0382016937182274          4.96996020133589         0.234703061000441
            -0.050835785331156          4.95235413212319          0.31171131955713
            -0.063368837475438          4.93107112088464         0.387802616532064
            -0.080031331450184          4.89677997700364         0.488496527009993
           -0.0964009322497484          4.85624298958237           0.5868917512885
            -0.112421694416715          4.80973573807352         0.682669598430554
            -0.128043214006999          4.75755396631082         0.775547365552795
            -0.145049134381423          4.69259930016282         0.876053292077056
            -0.161440622090589          4.62133675531964         0.972297821512827
            -0.177170855188286           4.5442183400676          1.06403463403466
            -0.192200683898424           4.4616921833756          1.15106418918647
            -0.208460001430162          4.36144550729978          1.24445379619915
            -0.223727676993895          4.25535741652712          1.33133420164216
            -0.237974755820245          4.14403024806315          1.41158207592263
            -0.251179996740741          4.02803600861358          1.48511868904983
            -0.264998508121562          3.89022745495227          1.56100635329738
            -0.277420062723431          3.74776784840078           1.6280359204465
            -0.288438334347593           3.6013626098573          1.68623106119387
             -0.29805302396734          3.45166893417197          1.73564872514799
            -0.305690157226073          3.31102513805948          1.77356512127723
            -0.312140099779861          3.16854541161365          1.80414205661417
            -0.317410841454917          3.02465155029153          1.82746590347732
            -0.321512191699292          2.87974094859598          1.84363335744803
            -0.324086010120158          2.75611838593223          1.85182325546023
             -0.32583103764415          2.63225039439751          1.85498666009876
            -0.326754775239335          2.50834604654139          1.85318935938568
            -0.326865142915572          2.38460586670073          1.84650016908841
            -0.326275533896468          2.27350659533658          1.83635200984885
            -0.325039226278059          2.16283268237184           1.8223458795187
            -0.323162267953922          2.05271649667228          1.80453239114339
            -0.320650904645271          1.94328656334394          1.78296380592703
            -0.317816162406802          1.84429695166026          1.76009026628178
            -0.314465046342121          1.74607420544357          1.73418168935613
            -0.310602737431102          1.64871039551544          1.70527578210538
            -0.306234841950284          1.55229604902619          1.67341092639009
            -0.301891488862627          1.46668910285076          1.64233987337137
            -0.297151277885037          1.38198403819887          1.60894200794167
            -0.292020488309853          1.29824754953426          1.57323919478166
            -0.286506942430475          1.21554731743889          1.53525146720517
             -0.28138366199999           1.1441871284886          1.50019992759481
            -0.275984641592807          1.07372234817198          1.46342060697659
            -0.270322637592699          1.00420905323808          1.42491126461638
            -0.264415150071804         0.935709522066566          1.38466207696518
            -0.260674671786117         0.894166217375383          1.35906776081747
            -0.256858008377776         0.853052131865475          1.33280039287953
            -0.252974834111282         0.812390245214114           1.3058479027292
            -0.249037352399907         0.772207291665762           1.2781943194599
            -0.245060896639183         0.732534752594493          1.24981892994302
            -0.241064920621637         0.693410683923897          1.22069440710064
            -0.237075071594782         0.654882889667731          1.19078498522456
            -0.233124151975677         0.617010363199244          1.16004563392081
            -0.230648106769647         0.593277684729897          1.14002213452332
            -0.228217811254025         0.569862973428515          1.11962216626016
            -0.225849195730495         0.546793725961579          1.09882530631869
            -0.223561134651479         0.524102360794245          1.07760816824531
            -0.221375979765588         0.501827182070164          1.05594408912478
            -0.219320384020607         0.480013900494922          1.03380236019565
             -0.21742722897067         0.458718962377532          1.01114776084701
            -0.215736525756438         0.438011116229606          0.98794034692489
            -0.214815006279457         0.425562280454738         0.973325545837781
            -0.214002732782736         0.413392305406461         0.958471302518458
            -0.213315821787941         0.401528961430753         0.943364905761404
            -0.212772997626917         0.390004577463667         0.927992816207752
            -0.212395871859393         0.378856684675315         0.912340849530155
            -0.212210123235207         0.368129816094555         0.896393895420993
            -0.212245380332493         0.357876560951221         0.880137586268896
            -0.212535068534041          0.34815795989262         0.863559130239709
            -0.212555269101556         0.347721432886043         0.862782455717491
            -0.212576092809974         0.347286203556321         0.862005068026316
            -0.212597544463318         0.346852280560303         0.861226966041292
            -0.212619628805482         0.346419672517159         0.860448148699362
            etc
            maglines{1} = 
                             0                         5                         0
          2.04176148777136e-05          5.00000675075849     -0.000126191739787496
          4.08351863333948e-05          5.00001349152292     -0.000252384020846009
          6.12527142408247e-05           5.0000202222932     -0.000378576842407675
          8.16701984737827e-05          5.00002694306922     -0.000504770203704611
          0.000183756960103661          5.00006039703174      -0.00113574507934097
          0.000285842610937721          5.00009360112219       -0.0017667333331647
          0.000387927135196574           5.0001265553277      -0.00239773486916731
           0.00049001051709984          5.00015925963549      -0.00302874959132762
           0.00100040973890324           5.0003190322655      -0.00618401763288604
           0.00151077803420038          5.00047255557266      -0.00933960091797448
           0.00202111342939052           5.0006198280352       -0.0124954874300385
           0.00253141395032944          5.00076084817451       -0.0156516651450471
            0.0050823243205154          5.00137211531144       -0.0314365004741103
             0.007632066440316          5.00182687747868       -0.0472268082115354
            0.0101803931388423          5.00212499876727        -0.063021078190313
            0.0127270571529473          5.00226637063118       -0.0788177970549258
            0.0254268141087839          5.00062008575691        -0.157785168816157
            0.0380479915269507          4.99505471881202        -0.236585705715436
             0.050560236463091          4.98558773348611         -0.31503067361234
            0.0629339077234446          4.97225354682649        -0.392934088443422
             0.080248929634465          4.94671071339667        -0.502515906401727
            0.0971470663849215          4.91364177005844        -0.610131749881702
             0.113556778167823          4.87330080668435        -0.715316060914159
             0.129413855219582          4.82598857790008        -0.817643835689873
             0.145166366572861          4.77010901873909        -0.920024974931304
             0.160214660805378          4.70754996455918           -1.018577661782
             0.174514423137449          4.63873590621319         -1.11298697911288
             0.188029865435701          4.56410188039581         -1.20298934934476
             0.201584892889142          4.47834162636552         -1.29412437291756
             0.214186420038259          4.38695848884691         -1.37976821150903
             0.225817342175233          4.29046844233697         -1.45976365693787
             0.236467471167393          4.18936941920011         -1.53399767246789
             0.247565765675621          4.06721322184781         -1.61263826361503
             0.257346238274611          3.94021914097906         -1.68340014317062
              0.26581681756111          3.80904279566006         -1.74626443525492
             0.272990831111069          3.67429746910476         -1.80124932365792
              0.27931986030199          3.52497111889617         -1.85197892153463
             0.284174467937203          3.37280010518301         -1.89360928244885
             0.287583415603357           3.2184105448381         -1.92624910758324
             0.289577965351715          3.06238521574017         -1.95002589158266
             0.290189970027169           2.9231623058127         -1.96380447895722
             0.289739631072324          2.78343426598698          -1.9708298572199
             0.288250342320926            2.643537038701         -1.97120464203313
             0.285746100721098          2.50378982173835         -1.96503561694684
             0.282782474631673          2.38335506859107         -1.95451447250768
             0.279093954618238          2.26344868065814         -1.93925251820126
             0.274695970532787          2.14425418627424         -1.91931726111897
             0.269604266565227          2.02594970844926         -1.89477741631806
             0.264419309738968          1.91994643953832         -1.86869669558625
             0.258692183506831          1.81493844330785         -1.83895912687854
             0.252434508134448          1.71105083902384         -1.80561255651647
             0.245658269138423          1.60840699412465         -1.76870596430945
              0.23912410711961          1.51716148547015          -1.7324741891927
             0.232188237973151          1.42711347059955         -1.69343023990886
             0.224859636629008          1.33835313568911         -1.65160855162918
             0.217147690103465          1.25097001296332         -1.60704552954098
             0.209977105527954          1.17450891280742         -1.56515431170683
             0.202517763296705          1.09927126652122         -1.52115162821083
             0.194776439370927           1.0253205326502         -1.47506558104857
             0.186760143415588         0.952718263190459         -1.42692814847088
              0.17978738797092         0.892550938237161           -1.384741413986
             0.172627591827363         0.833421496160737         -1.34114233755801
             0.165283740846495         0.775361756433708         -1.29615833193343
              0.15775769951305         0.718398386415769         -1.24982287550136
             0.153294691396793         0.685766162758762         -1.22225497535831
             0.148769160676012         0.653515463727945         -1.19425174542332
             0.144179283610778         0.621647981103645         -1.16582348908377
             0.139521848972148         0.590163927723201         -1.13698196036292
             0.135852118532761         0.565957860601514         -1.11428863286708
             0.132134836427589         0.541981425888754         -1.09136047578431
             0.128364354097521         0.518232891519742         -1.06820494155939
             0.124532555291891         0.494710113784141         -1.04483006830521
             0.121496076465136         0.476529511665792         -1.02646300922512
             0.118408229959158         0.458483701280984         -1.00797198952587
             0.115259030554191         0.440572207038815        -0.989361112660354
             0.112035290043433         0.422795273698917        -0.970634359085049
             0.108719770378775         0.405154287736759         -0.95179545237296
              0.10529020377374         0.387652246897418        -0.932847598383982
             0.101716030257482         0.370295456055628        -0.913793186729339
            0.0979567899708584         0.353094357349851         -0.89463364431595
            0.0972485225864908         0.349980088145269        -0.891133641177215
            0.0965320260395926         0.346871617282926        -0.887630165258499
            0.0958069208751358         0.343769067108972        -0.884123211821123
            0.0950727966258886         0.340672575791845          -0.8806127755025
            etc
