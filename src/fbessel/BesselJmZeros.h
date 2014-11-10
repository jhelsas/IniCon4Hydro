#define MaxNumOfBesselJmIndex 21
#define MaxNumOfBesselJmZeros 50
double lambda[MaxNumOfBesselJmIndex][MaxNumOfBesselJmZeros] = {
  {2.404825557695770,5.520078110286310,8.653727912911010,11.791534439014300,14.930917708487800,18.071063967910900,21.211636629879300,24.352471530749300,27.493479132040200,30.634606468432000,33.775820213573600,36.917098353664000,40.058425764628200,43.199791713176700,46.341188371661800,49.482609897397800,52.624051841115000,55.765510755020000,58.906983926080900,62.048469190227200,65.189964800206900,68.331469329856800,71.472981603593700,74.614500643701900,77.756025630388100,80.897555871137600,84.039090776938200,87.180629843641200,90.322172637210500,93.463718781944800,96.605267950996300,99.746819858680600,102.888374254195000,106.029930916452000,109.171489649805000,112.313050280495000,115.454612653667000,118.596176630873000,121.737742087951000,124.879308913233000,128.020877006008000,131.162446275214000,134.304016638306000,137.445588020284000,140.587160352854000,143.728733573690000,146.870307625797000,150.011882456955000,153.153458019228000,156.295034268534000},
  {3.831705970207510,7.015586669815610,10.173468135062700,13.323691936314200,16.470630050877600,19.615858510468200,22.760084380592800,25.903672087618400,29.046828534916900,32.189679910974400,35.332307550083900,38.474766234771600,41.617094212814500,44.759318997652800,47.901460887185500,51.043535183571500,54.185553641061300,57.327525437901000,60.469457845347500,63.611356698481200,66.753226734098500,69.895071837495800,73.036895225573800,76.178699584641500,79.320487175476300,82.462259914373600,85.604019436350200,88.745767144926300,91.887504251695000,95.029231808044700,98.170950730790800,101.312661823039000,104.454365791283000,107.596063259509000,110.737754780899000,113.879440847595000,117.021121898892000,120.162798328149000,123.304470488636000,126.446138698517000,129.587803245104000,132.729464388510000,135.871122364789000,139.012777388660000,142.154429655859000,145.296079345196000,148.437726620342000,151.579371631401000,154.721014516286000,157.862655401930000},
  {5.135622301840680,8.417244140399860,11.619841172149100,14.795951782351300,17.959819494987800,21.116997053021900,24.270112313573100,27.420573549984600,30.569204495516400,33.716519509222700,36.862856511283800,40.008446733478200,43.153453778371500,46.297996677236900,49.442164110416900,52.586023506816000,55.729627053201100,58.873015772612200,62.016222359217700,65.159273190757800,68.302189784183500,71.444989866357800,74.587688173602400,77.730297056978900,80.872826946244800,84.015286709546200,87.157683935203400,90.300025154592900,93.442316020011100,96.584561447783200,99.726765734292800,102.868932650728000,106.011065520963000,109.153167285982000,112.295240557472000,115.437287662664000,118.579310682042000,121.721311481196000,124.863291737881000,128.005252965073000,131.147196530718000,134.289123674703000,137.431035523503000,140.572933102855000,143.714817348778000,146.856689117169000,149.998549192200000,153.140398293676000,156.282237083508000,159.424066171418000},
  {6.380161895923980,9.761023129981670,13.015200721698400,16.223466160318800,19.409415226435000,22.582729593104400,25.748166699295000,28.908350780921800,32.064852407097700,35.218670738610100,38.370472434756900,41.520719670406800,44.669743116617200,47.817785691533300,50.965029906205200,54.111615569821900,57.257651604499000,60.403224138472100,63.548402178567200,66.693241667372700,69.837788437904300,72.982080400432000,76.126149184774100,79.270021390055900,82.413719547267900,85.557262868830000,88.700667838222100,91.843948678147100,94.987117725465600,98.130185733874900,101.273162120080000,104.416055165397000,107.558872181933000,110.701619650395000,113.844303335032000,116.986928380009000,120.129499390634000,123.272020502131000,126.414495438148000,129.556927560730000,132.699319913184000,135.841675256988000,138.983996103681000,142.126284742523000,145.268543264558000,148.410773583617000,151.552977454706000,154.695156490148000,157.837312173800000,160.979445873597000},
  {7.588342434503800,11.064709488501200,14.372536671617600,17.615966049804800,20.826932956962400,24.019019524771100,27.199087765981300,30.371007667117200,33.537137711819200,36.699001128744600,39.857627302180900,43.013737723354400,46.167853512924400,49.320360686390300,52.471551398458000,55.621650909768000,58.770835740459200,61.919246204097700,65.066995255695600,68.214174861467000,71.360860665298000,74.507115461396400,77.652991815343400,80.798534067923700,83.943779885098100,87.088761469813600,90.233506518792400,93.378038984848900,96.522379689381200,99.666546818328800,102.810556326690000,105.954422270697000,109.098157082311000,112.241771797409000,115.385276246530000,118.528679215153000,121.671988579057000,124.815211419159000,127.958354119370000,131.101422450345000,134.244421641424000,137.387356442681000,140.530231178623000,143.673049794826000,146.815815898563000,149.958532794316000,153.101203514908000,156.243830848863000,159.386417364526000,162.528965431386000},
  {8.771483815959950,12.338604197466900,15.700174079711700,18.980133875179900,22.217799896561300,25.430341154222700,28.626618307291100,31.811716724047800,34.988781294559300,38.159868561967100,41.326383254047400,44.489319123219700,47.649399806697000,50.807165203006300,53.963026558378100,57.117302781504200,60.270245072942800,63.422054045875800,66.572891887118300,69.722891161716800,72.872161296912000,76.020793430591600,79.168864087087400,82.316437999012300,85.463570298373100,88.610308235796300,91.756692542506100,94.902758518889600,98.048536911695200,101.194054626309000,104.339335309226000,107.484399827534000,110.629266666070000,113.773952258293000,116.918471263455000,120.062836799964000,123.207060642831000,126.351153391488000,129.495124613039000,132.638982965051000,135.782736301208000,138.926391762565000,142.069955856637000,145.213434526178000,148.356833209185000,151.500156891406000,154.643410152421000,157.786597206202000,160.929721936904000,164.072787930528000},
  {9.936109524217680,13.589290170541200,17.003819667816000,20.320789213566500,23.586084435581400,26.820151983411400,30.033722386570500,33.233041762847100,36.422019668258500,39.603239416075400,42.778481613199500,45.949015998042600,49.115773724764300,52.279453903601100,55.440592068853200,58.599605631237700,61.756824901876800,64.912514784720700,68.066890268038700,71.220127696168400,74.372373108624300,77.523748502423500,80.674356598679300,83.824284515391900,86.973606629195900,90.122386828076200,93.270680301413900,96.418534974778200,99.565992669243900,102.713090045137000,105.859859375645000,109.006329185097000,112.152524778781000,115.298468685253000,118.444181027558000,121.589679836376000,124.734981315426000,127.880100067424000,131.025049287279000,134.169840927934000,137.314485843293000,140.458993911817000,143.603374143806000,146.747634774809000,149.891783347210000,153.035826781723000,156.179771440193000,159.323623180944000,162.467387407659000,165.611069112680000},
  {11.086370019245100,14.821268727013200,18.287582832481700,21.641541019848400,24.934927887673000,28.191188459483200,31.422794192265600,34.637089352069300,37.838717382853600,41.030773691585500,44.215408505261300,47.394165755570500,50.568184679795600,53.738325371963300,56.905249991978800,60.069476998277000,63.231418368888300,66.391405759533000,69.549709272422300,72.706551172477200,75.862116076322400,79.016558632922400,82.170009390527900,85.322579332379500,88.474363421866700,91.625443401413900,94.775890022676600,97.925764838797900,101.075121656129000,104.224007718760000,107.372464681664000,110.520529415269000,113.668234674663000,116.815609659307000,119.962680483652000,123.109470574803000,126.256001010118000,129.402290805071000,132.548357159757000,135.694215670794000,138.839880514185000,141.985364603684000,145.130679728403000,148.275836672787000,151.420845321529000,154.565714751593000,157.710453313159000,160.855068701009000,163.999568017662000,167.143957829332000},
  {12.225092264004600,16.037774190887700,19.554536430997100,22.945173131874600,26.266814641176600,29.545659670998500,32.795800037341500,36.025615063869600,39.240447995178100,42.443887743273600,45.638444182199100,48.825930381553900,52.007691456686900,55.184747939289100,58.357889025269700,61.527735166816000,64.694781235818700,67.859426993000800,71.021999040620500,74.182766927652800,77.341955156796000,80.499752266331700,83.656317789561700,86.811787651267100,89.966278397575300,93.119890544332300,96.272711251866000,99.424816479638900,102.576272735453000,105.727138505778000,108.877465433208000,112.027299291864000,115.176680800158000,118.325646301789000,121.474228339291000,124.622456139452000,127.770356026037000,130.917951772231000,134.065264902867000,137.212314954575000,140.359119700577000,143.505695345599000,146.652056695447000,149.798217305015000,152.944189607862000,156.089985029986000,159.235614090000000,162.381086487578000,165.526411181733000,168.671596460278000},
  {13.354300477435200,17.241220382489100,20.807047789264100,24.233885257750600,27.583748963573000,30.885378967696700,34.154377923855100,37.400099977156600,40.628553718964500,43.843801420337400,47.048700737654000,50.245326955305400,53.435227157042100,56.619580266508400,59.799301630960200,62.975113534241600,66.147594024796000,69.317211517895100,72.484349817473000,75.649326536060800,78.812406871964200,81.973814061805500,85.133737413339200,88.292338551132200,91.449756324634400,94.606110702857700,97.761505892707800,100.916032856444000,104.069771359662000,107.222791649243000,110.375155837254000,113.526919049415000,116.678130383716000,119.828833714929000,122.979068373244000,126.128869719480000,129.278269634861000,132.427296939834000,135.575977753672000,138.724335804423000,141.872392697048000,145.020168146185000,148.167680178897000,151.314945311807000,154.461978706356000,157.608794305243000,160.755404952694000,163.901822500729000,167.048057903314000,170.194121299969000},
  {14.475500686554700,18.433463666966600,22.046985364697800,25.509450554182800,28.887375063530500,32.211856199712700,35.499909205373800,38.761807017881700,42.004190236671800,45.231574103535000,48.447151387269400,51.653251668165900,54.851619075963400,58.043587928232500,61.230197977292700,64.412272412924400,67.590472073698500,70.765333996242800,73.937299381768100,77.106734246861300,80.273944913985200,83.439189796105700,86.602688476727600,89.764628787179100,92.925172381168400,96.084459168145400,99.242610870419200,102.399733900616000,105.555921706997000,108.711256698525000,111.865811835334000,115.019651950819000,118.172834856929000,121.325412273188000,124.477430611490000,127.628931642237000,130.779953062279000,133.930528981211000,137.080690339417000,140.230465268832000,143.379879405383000,146.528956160520000,149.677716957976000,152.826181440835000,155.974367653194000,159.122292199979000,162.269970387923000,165.417416350266000,168.564643157316000,171.711662914721000},
  {15.589847884455500,19.615966903966900,23.275853726263400,26.773322545509500,30.179061178784900,33.526364075588600,36.833571341894900,40.111823270954200,43.368360947521700,46.608132676274900,49.834653510396700,53.050498959135000,56.257604715114500,59.457456908388000,62.651217388202900,65.839808804445100,69.023973928461300,72.204317964118300,75.381339332109500,78.555452463758800,81.727004943537600,84.896290582739800,88.063559516403000,91.229026090787500,94.392875089319200,97.555266694032800,100.716340474093000,103.876218618181000,107.035008573751000,110.192805217004000,113.349692648597000,116.505745688641000,119.661031128389000,122.815608783802000,125.969532386768000,129.122850342554000,132.275606376414000,135.427840087909000,138.579587427992000,141.730881111163000,144.881750972819000,148.032224280137000,151.182326003401000,154.332079053556000,157.481504490790000,160.630621708210000,163.779448594004000,166.928001675003000,170.076296244067000,173.224346473403000},
  {16.698249933848200,20.789906360078400,24.494885043881300,28.026709949973100,31.459960035318000,34.829986990290200,38.156377504681400,41.451092307939700,44.721943543191200,47.974293531269100,51.211967004101100,54.437776928325100,57.653844811907000,60.861804682480500,64.062937824850100,67.258264556341300,70.448608396520600,73.634641960749100,76.816920437231600,79.995906435637300,83.171988718832900,86.345496520534100,89.516710626564300,92.685872048911200,95.853188885871900,99.018841799052700,102.182988424175000,105.345766951770000,108.507299055586000,111.667692304064000,114.827042158870000,117.985433641143000,121.142942728486000,124.299637532387000,127.455579295476000,130.610823240133000,133.765419293776000,136.919412711324000,140.072844611519000,143.225752440751000,146.378170375608000,149.530129673427000,152.681658978539000,155.832784590626000,158.983530700575000,162.133919598331000,165.283971856569000,168.433706493424000,171.583141116995000,174.732292053997000},
  {17.801435153282400,21.956244067836300,25.705103053924700,29.270630441874800,32.731053310978400,36.123657666448800,39.469206825243900,42.780439265447200,46.065710911575600,49.330780096443500,52.579769064383400,55.815719876305800,59.040934037249300,62.257189393731700,65.465883797232100,68.668133216891100,71.864840513502700,75.056744744974100,78.244457215495800,81.428488292663100,84.609267666025000,87.787159863201400,90.962476282045300,94.135484626561100,97.306416382892700,100.475472798191000,103.642829703475000,106.808641434976000,109.973044045975000,113.136157955500000,116.298090146509000,119.458936001046000,122.618780840814000,125.777701227222000,128.935766063817000,132.093037535460000,135.249571911910000,138.405420238223000,141.560628930212000,144.715240289938000,147.869292953519000,151.022822281465000,154.175860699977000,157.328438000290000,160.480581601968000,163.632316785142000,166.783666895884000,169.934653528305000,173.085296686382000,176.235614928125000},
  {18.899997953174000,23.115778347252800,26.907368976182100,30.505950163896000,33.993184984781500,37.408185128639700,40.772827853501900,44.100590565798300,47.400347780543200,50.678236946479900,53.938666209126900,57.184898598119300,60.419409852130300,63.644117508962300,66.860533012260100,70.069865833196500,73.273096621009900,76.471029759814700,79.664331875077600,82.853560536123500,86.039185980685900,89.221607784583600,92.401167811358000,95.578160384997100,98.752840362876500,101.925429602183000,105.096122183942000,108.265088666711000,111.432479575495000,114.598428282808000,117.763053402792000,120.926460792437000,124.088745233582000,127.249991853948000,130.410277333505000,133.569670933314000,136.728235376726000,139.886027607239000,143.043099442769000,146.199498142580000,149.355266900237000,152.510445273653000,155.665069561443000,158.819173133279000,161.972786720694000,165.125938673781000,168.278655188379000,171.430960507638000,174.582877101291000,177.734425825467000},
  {19.994430629816400,24.269180026208900,28.102415231667800,31.733413344374600,35.247086785793300,38.684276389289600,42.067916998656800,45.412189614733100,48.726464116240700,52.017241278881600,55.289204146560000,58.545828904385100,61.789759895945000,65.023050251042200,68.247321996420800,71.463875885022700,74.673768712140400,77.877868973486400,81.076897720632800,84.271459069716400,87.462063333480800,90.649144800867800,93.833075571297500,97.014176439321700,100.192725545543000,103.368965316046000,106.543108076381000,109.715340628929000,112.885828012170000,116.054716608866000,119.222136732042000,122.388204789077000,125.553025102660000,128.716691450866000,131.879288375967000,135.040892301757000,138.201572491493000,141.361391872508000,144.520407748781000,147.678672418929000,150.836233714011000,153.993135467106000,157.149417924579000,160.305118107373000,163.460270129285000,166.614905478118000,169.769053264692000,172.922740443919000,176.075992011565000,179.228831179766000},
  {21.085146113064700,25.417019006342800,29.290870696313400,32.953664885068700,36.493397912446500,39.952553490221200,43.355073203851900,46.715809435025800,50.044606016842500,53.348312331852900,56.631875942813000,59.898978728776800,63.152428193380300,66.394409038588300,69.626650879869700,72.850543506754200,76.067218169773300,79.277606200069500,82.482482112139900,85.682495843831100,88.878197240069300,92.070054900383000,95.258470865972800,98.443792191938400,101.626320157431000,104.806317663574000,107.984015226107000,111.159615867587000,114.333299140072000,117.505224454971000,120.675533856560000,123.844354345550000,127.011799836273000,130.177972813670000,133.342965742841000,136.506862273506000,139.669738273616000,142.831662719890000,145.992698468022000,149.152902921214000,152.312328612441000,155.471023713244000,158.629032479689000,161.786395644421000,164.943150762292000,168.099332515884000,171.254972986293000,174.410101893692000,177.564746811572000,180.718933357975000},
  {22.172494618826300,26.559784138025400,30.473279946333500,34.167267853840500,37.732680522054100,41.213567059135000,44.634829753112700,48.011962936065500,51.355264650751700,54.671919177579500,57.967128833466200,61.244774098170400,64.507820402735200,67.758580113886900,70.998887489854200,74.230219119065500,77.453778998172600,80.670559983454700,83.881389045812700,87.086961171532000,90.287865144784200,93.484603422990700,96.677607646035700,99.867250872460300,103.053857330537000,106.237710260391000,109.419058274081000,112.598120553797000,115.775091130951000,118.950142432118000,122.123428235648000,125.295086151153000,128.465239710077000,131.634000137300000,134.801467859588000,137.967733795719000,141.132880464556000,144.296982940563000,147.460109680874000,150.622323243749000,153.783680914811000,156.944235254662000,160.104034579227000,163.263123382324000,166.421542708460000,169.579330482586000,172.736521802536000,175.893149199013000,179.049242867252000,182.204830873932000},
  {23.256776085110000,27.697898350855300,31.650118151857000,35.374717217954900,38.965432047654000,42.467807213330800,45.907663866365000,49.301111338030600,52.658883651520200,55.988487220554300,59.295369944286600,62.583604180135600,65.856308280514900,69.115918502286000,72.364370871494500,75.603226568097300,78.833760629364900,82.057026112746600,85.273901414437400,88.485125764711100,91.691326259226000,94.893038724129100,98.090724018109100,101.284780909837000,104.475556352223000,107.663353754797000,110.848439700215000,114.031049439705000,117.211391421571000,120.389651047640000,123.565993808475000,126.740567915136000,129.913506520208000,133.084929601646000,136.254945568175000,139.423652633513000,142.591139997636000,145.757488866209000,148.922773333676000,152.087061150945000,155.250414395030000,158.412890055035000,161.574540546503000,164.735414164215000,167.895555481901000,171.055005706040000,174.213802989810000,177.371982712360000,180.529577727821000,183.686618587818000},
  {24.338249623407200,28.831730351300900,32.821802761873600,36.576450758999600,40.192095100837100,43.715712417955400,47.174004568236500,50.583671140024000,53.955865282862700,57.298403654373200,60.616971127173500,63.915825576156300,67.198233500661600,70.466751417355000,73.723414330835800,76.969865851345900,80.207450370705300,83.437279825079500,86.660282994477400,89.877242530642300,93.088823190097400,96.295593652912800,99.498043589692800,102.696597158509000,105.891623785257000,109.083446852794000,112.272350763137000,115.458586721572000,118.642377507701000,121.823921436827000,125.003395669271000,128.180958990760000,131.356754160934000,134.530909906992000,137.703542624063000,140.874757831883000,144.044651427902000,147.213310769536000,150.380815612359000,153.547238926266000,156.712647607903000,159.877103104502000,163.040661961827000,166.203376306842000,169.365294274065000,172.526460383176000,175.686915874286000,178.846699006352000,182.005845323398000,185.164387892543000},
  {25.417140814072500,29.961603791625200,33.988702785235200,37.772857844399000,41.413065513892600,44.957676748421700,48.434239195205700,51.860019928074600,55.246575614639300,58.602022073846700,61.932273072882600,65.241765993816600,68.533910938821000,71.811381203719600,75.076308077035800,78.330415494843200,81.575115548130200,84.811577742024900,88.040780199243000,91.263548162504400,94.480583385064500,97.692486869060200,100.899776670940000,104.102901997134000,107.302254474131000,110.498177241360000,113.690972348690000,116.880906820843000,120.068217664244000,123.253116027873000,126.435790682236000,129.616410944746000,132.795129152718000,135.972082764368000,139.147396152136000,142.321182140135000,145.493543327729000,148.664573233451000,151.834357287339000,155.002973694796000,158.170494191134000,161.336984702723000,164.502505928064000,167.667113849956000,170.830860188168000,173.993792800593000,177.155956039619000,180.317391069506000,183.478136149669000,186.638226888093000}
};