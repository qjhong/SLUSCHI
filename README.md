<!-- wp:wp-bootstrap-blocks/container -->
<!-- wp:wp-bootstrap-blocks/row {"template":"custom"} -->
<!-- wp:wp-bootstrap-blocks/column {"sizeMd":10} -->


<!-- wp:heading -->
<h2 class="wp-block-heading">Documentation</h2>
<!-- /wp:heading -->

<!-- wp:paragraph -->
<p><a href="http://creativecommons.org/licenses/by-nd/4.0/"></a><br>This work is licensed under a <a href="http://creativecommons.org/licenses/by-nd/4.0/">Creative Commons Attribution-NoDerivatives 4.0 International License</a>.</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>Please contact me at qhong@alumni.caltech.edu if you have questions.</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>Several Notes for Users:</p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true} -->
<ol class="uds-list"><!-- wp:list-item -->
<li>A typical calculation with SLUSCHI takes a few days/weeks, with CPU cost ranging from 5,000 to 100,000 CPU hours. On <a href="[https://blogs.brown.edu/qhong/?page_id=7](https://faculty.engineering.asu.edu/hong/systems-tested/)">this page</a>, I list the systems I have studied.</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Try "kmesh = -1" in "job.in". This activates the use of a special kpoint (1/4,1/4,1/4). This provides a low-cost option (only twice the cost of gamma version), and yet it often generates reliable results for melting point.</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>An example of SLUSCHI run is available at the <a href="https://materialsdata.nist.gov/handle/11256/693">repositories</a> at the National Institute of Standards and Technology. Or follow this <a href="https://drive.google.com/drive/folders/1DuiBRLLoJCH1EfLM0qCZyJsuIUR62S6l?usp=sharing">link</a> on Google Drive for the same copy.</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
<p><p>Sluschi_2.0 <br>What's new:</p></p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true} -->
<ol class="uds-list"><!-- wp:list-item -->
<li>mds: new capability to compute entropy of both solids and liquids from one single MD trajectory. User manual coming soon. </li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Bug fixes and various improvements.</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
<p><p>Sluschi_1.3<br>What's new:</p></p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true} -->
<ol class="uds-list"><!-- wp:list-item -->
<li>Add feature to detect failed VASP jobs and restart them. It takes hundreds of VASP MD runs to calculate a melting temperature, and some of them may accidentally fail (e.g., running beyond walltime limit, failed to start due to queue/disk/network issue, etc.). Now SLUSCHI checks VASP MD runs and restarts failed jobs. Set "detectfail" and "maxwaithour" in job.in.</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Bug fixes and various improvements.</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
<p>Sluschi_1.2<br>What's new:</p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true} -->
<ol class="uds-list"><!-- wp:list-item -->
<li>Bug fixes and various improvements.</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
<p>Sluschi_1.1<br>What's new:</p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true} -->
<ol class="uds-list"><!-- wp:list-item -->
<li>"Bold sampling" mode: SLUSCHI launches MD duplicates aggressively. This strategy significantly reduces the physical time of calculations, though it may slightly increase computer cost.</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Heat of fusion calculations</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Bug fixes and various improvements.</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:paragraph -->
<p><a href="https://brownbox.brown.edu/download.php?hash=b9fff340">Sluschi_1</a><br>User guide at <a href="https://www.sciencedirect.com/science/article/pii/S0364591615300468">the <em>CALPHAD Journal</em></a>.<br>What's new:</p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true} -->
<ol class="uds-list"><!-- wp:list-item -->
<li>intpol: user does not need to guess a melting point. See "intpol".</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>adj_bmix: let SLUSCHI decide the value of BMIX.</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Bug fixes and various improvements.</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:paragraph {"className":"lead mt-lg-12 mt-6"} -->
<p class="lead mt-lg-12 mt-6">SLUSCHI is a fully automated code which calculates melting points based on first-principles molecular dynamics simulations, with interface to the first-principles code VASP. Starting from the crystal structure of a solid (which the user inputs), SLUSCHI will automatically build a supercell of a proper size, prepare solid-liquid coexistence, and then employ the small-cell coexistence method to calculate the melting temperature. SLUSCHI is applicable to a wide variety of materials, thanks to the fact that density functional theory calculations are highly generalizable.</p>
<!-- /wp:paragraph -->

<!-- wp:image {"id":83,"sizeSlug":"large","linkDestination":"none"} -->
<figure class="wp-block-image size-large"><img src="https://faculty.engineering.asu.edu/hong/wp-content/uploads/sites/65/2015/02/Methods-1500x1369.jpg" alt="" class="wp-image-83"/></figure>
<!-- /wp:image -->

<!-- wp:separator {"opacity":"css","className":"mt-12 is-style-copy-divider"} -->
<hr class="wp-block-separator has-css-opacity mt-12 is-style-copy-divider"/>
<!-- /wp:separator -->

<!-- wp:heading -->
<h2 class="wp-block-heading">Systems Tested</h2>
<!-- /wp:heading -->

<!-- wp:paragraph -->
<p>Ta (0 and 200 GPa, bcc), Na (30-120 GPa, bcc/fcc), NaCl, La<sub>2</sub>Zr<sub>2</sub>O<sub>7</sub> (La<sub>2</sub>O<sub>3</sub>-2ZrO<sub>2</sub>), the Hf-Ta-C-N systems, Al, Si, HfO<sub>2</sub>, Fe (330 GPa), Ti (bcc) and many more…</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>m.p.: melting point calculated by SLUSCHI.<br>If not specified, PAW-PBE is employed in the calculations.<br>HSE: melting point after HSE correction. (important for semiconductors)</p>
<!-- /wp:paragraph -->

<!-- wp:heading {"level":4} -->
<h4 class="wp-block-heading">MD: total number of MD trajectories</h4>
<!-- /wp:heading -->

<!-- wp:paragraph -->
<p>CPU hours are measured as / converted to TACC Stampede.<br>days: the physical time it takes. SLUSCHI is currently optimized to lower the CPU cost, rather the physical time. In order to reduce the physical time, user may manually explore temperature of interest. However, this may increase CPU hours.<br>"early runs" are calculations performed at early stage of method development. Hence the efficiency is relatively low compared to "sluschi".</p>
<!-- /wp:paragraph -->

<!-- wp:paragraph -->
<p>poor results and reason&nbsp;|&nbsp;good results</p>
<!-- /wp:paragraph -->

<!-- wp:table {"className":"is-style-stripes"} -->
<figure class="wp-block-table is-style-stripes"><table><tbody><tr><td>systems</td><td>m.p./K</td><td>HSE/K</td><td>expt./K</td><td>DFT&nbsp;PP</td><td>rad / Å</td><td>kmesh</td><td># MD</td><td>cpu hours</td><td>days</td><td>note</td></tr><tr><td>Al</td><td>1040±13</td><td></td><td>933</td><td>Al</td><td>12</td><td>(1/2,1/2,1/2)</td><td>19</td><td>5,400</td><td>7</td><td>sluschi</td></tr><tr><td>Al</td><td>999±21</td><td>&nbsp;1054</td><td>933</td><td>Al</td><td>10</td><td>Auto 30</td><td>58</td><td>16,000</td><td>15</td><td>sluschi</td></tr><tr><td>Ti_v</td><td>1811±47</td><td></td><td>1941</td><td>Ti_v</td><td>9</td><td>Auto 20</td><td>49</td><td>&nbsp;15,500</td><td>21</td><td>sluschi</td></tr><tr><td>Ti_v</td><td>1750±25</td><td>1971</td><td>1941</td><td>Ti_v</td><td>10</td><td>Auto 20</td><td>26</td><td>7,900</td><td>17</td><td>sluschi</td></tr><tr><td>Ti_pv</td><td>1952±45</td><td></td><td>1941</td><td>Ti_pv</td><td>10</td><td>(1/4,1/4,1/4)</td><td>&nbsp;48</td><td>20,000</td><td>20</td><td>sluschi</td></tr><tr><td>Ti_pv 60GPa</td><td>2505±47</td><td></td><td>?</td><td>Ti_pv</td><td>10</td><td>(1/4,1/4,1/4)</td><td>&nbsp;46</td><td>39,500</td><td>17</td><td>sluschi</td></tr><tr><td>Si</td><td>724±47</td><td></td><td>1687</td><td>Si</td><td>10</td><td>Auto 10</td><td>25</td><td>1,000</td><td>2</td><td>sluschi</td></tr><tr><td>Si</td><td>1378±24</td><td>1785</td><td>1687</td><td>Si</td><td>10</td><td>Auto 20</td><td>19</td><td>7,500</td><td>15</td><td>el. DOS change, require HSE</td></tr><tr><td>Si</td><td>1364±37</td><td></td><td>1687</td><td>Si</td><td>12</td><td>Auto 20</td><td>20</td><td>27,700</td><td>27</td><td>sluschi</td></tr><tr><td>diamond,100GPa</td><td>4307±8</td><td></td><td>4250</td><td>C</td><td>10</td><td>&nbsp;gamma</td><td>45</td><td>313,000</td><td>168</td><td>sluschi</td></tr><tr><td>Ru_v</td><td>2435±32</td><td></td><td>2607</td><td>Ru_v</td><td>10</td><td>&nbsp;Auto 20</td><td>&nbsp;30</td><td>88,000</td><td>37</td><td>sluschi</td></tr><tr><td>Ru_pv</td><td>2550±34</td><td></td><td>2607</td><td>Ru_pv</td><td>10</td><td>(1/4,1/4,1/4)</td><td>&nbsp;33</td><td>139,000</td><td>51</td><td>sluschi</td></tr><tr><td>Ru ternary alloys</td><td>2564±40</td><td></td><td>n/a</td><td>&nbsp;pv</td><td>10</td><td>(1/4,1/4,1/4)</td><td>&nbsp;23</td><td>247,000</td><td>80</td><td>sluschi</td></tr><tr><td>Hf,bcc</td><td>&nbsp;2562±31</td><td></td><td>2506</td><td>Hf</td><td>10</td><td>(1/2,1/2,1/2)</td><td>115</td><td>14,900</td><td>12</td><td>sluschi</td></tr><tr><td>Hf,hcp</td><td>&nbsp;2122±50</td><td></td><td>n/a</td><td>Hf</td><td>10</td><td>(1/2,1/2,1/2)</td><td>50</td><td>25,800</td><td>21</td><td>sluschi,&nbsp;hcp not stable at MT</td></tr><tr><td>HfO<sub>2</sub></td><td>&nbsp;2327±47</td><td>-</td><td>3031</td><td>valence</td><td>10</td><td>gamma</td><td>33</td><td>79,000</td><td></td><td>sluschi</td></tr><tr><td>HfO<sub>2</sub> PBE+U</td><td>&nbsp;3486±91</td><td>-</td><td>3031</td><td>Hf_pv</td><td>10</td><td>gamma</td><td>34</td><td>99,600</td><td></td><td>sluschi</td></tr><tr><td>HfO<sub>2</sub> PBE+U</td><td>&nbsp;3313±81</td><td>-</td><td>3031</td><td>Hf_pv</td><td>12</td><td>gamma</td><td>25</td><td>running</td><td></td><td>ionic, use&nbsp;large $rad</td></tr><tr><td>ZrO<sub>2</sub></td><td></td><td></td><td>2988</td><td>valence</td><td>10</td><td>Auto 20</td><td></td><td></td><td></td><td>running</td></tr><tr><td>Ta_v</td><td>2986±41</td><td></td><td>3290</td><td>Ta_v</td><td>10</td><td>(0,1/4,1/4)</td><td>38</td><td>32,000</td><td>23</td><td>early runs, low efficiency</td></tr><tr><td>Ta_pv</td><td>3194±40</td><td></td><td>3290</td><td>Ta_pv</td><td>10</td><td>(0,1/4,1/4)</td><td>38</td><td>54,000</td><td>24</td><td>early runs, low efficiency</td></tr><tr><td>Ta_pv_PW91</td><td>3066±51</td><td></td><td>3290</td><td>Ta_pv</td><td>10</td><td>(0,1/4,1/4)</td><td>38</td><td>38,000</td><td>68</td><td>PW91 [1]</td></tr><tr><td>Ta_pv,200GPa</td><td>7953±69</td><td></td><td>n/a</td><td>Ta_pv</td><td>10</td><td>(1/4,1/4,1/4)</td><td>80</td><td>150,000</td><td>48</td><td>sluschi, high efficiency [4]</td></tr><tr><td>W</td><td>3497±54</td><td></td><td>3695</td><td>W</td><td>10</td><td>A20</td><td>22</td><td>35,900</td><td>49</td><td>sluschi</td></tr><tr><td>W_pv</td><td>3470±45</td><td></td><td>3695</td><td>W_pv</td><td>10</td><td>(1/4,1/4,1/4)</td><td>30</td><td>38,500</td><td>18</td><td>sluschi</td></tr><tr><td>Na 15 GPa</td><td>657±8</td><td></td><td>810, 698 ?</td><td>Na_pv</td><td>10.4</td><td>(0,1/4,1/4)</td><td>55</td><td>47,000</td><td>24</td><td>bcc, expt under dispute, e.g.,</td></tr><tr><td>Na 26 GPa</td><td>750±16</td><td>706</td><td>970 ?</td><td>Na_pv</td><td>9.8</td><td>(0,1/4,1/4)</td><td>52</td><td>26,000</td><td>&nbsp;17</td><td>Zha, Boehler vs. Gregoryanz</td></tr><tr><td>Na 40 GPa</td><td>742±17</td><td></td><td>950 ?</td><td>Na_pv</td><td>9.4</td><td>(0,1/4,1/4)</td><td>74</td><td>54,000</td><td>37</td><td>SLUSCHI results&nbsp;agree with</td></tr><tr><td>Na 55 GPa</td><td>716±12</td><td></td><td>810 ?</td><td>Na_pv</td><td>9.0</td><td>(0,1/4,1/4)</td><td>56</td><td>64,000</td><td>&nbsp;31</td><td>Eshet and&nbsp;Desjarlais (theory).</td></tr><tr><td>Na 80 GPa</td><td>674±20</td><td></td><td>700 ?</td><td>Na_pv</td><td>10.9</td><td>(0,1/4,1/4)</td><td>55</td><td>170,000</td><td>67</td><td>fcc, expt under dispute</td></tr><tr><td>Na 100 GPa</td><td>662±14</td><td></td><td>450 ?</td><td>Na_pv</td><td>10.6</td><td>(0,1/4,1/4)</td><td>48</td><td>57,000</td><td>28</td><td>fcc, expt under dispute</td></tr><tr><td>Na 120 GPa</td><td>579±27</td><td></td><td>?</td><td>Na_pv</td><td>10.4</td><td>(0,1/4,1/4)</td><td>88</td><td>150,000</td><td>77</td><td>fcc, expt under dispute</td></tr><tr><td>NaCl</td><td>1014±18</td><td></td><td>1074</td><td>valence</td><td>11</td><td>gamma</td><td>59</td><td>24,000</td><td>50</td><td>early runs, low efficiency</td></tr><tr><td>La<sub>2</sub>Zr<sub>2</sub>O<sub>7</sub></td><td>2420±27</td><td>2630</td><td>2530</td><td>Zr_v, La_sv</td><td>10</td><td>gamma</td><td>64</td><td>&nbsp;575,000</td><td>210</td><td>[2]</td></tr><tr><td>Hf-Ta-C-N</td><td>-</td><td></td><td>-</td><td>valence</td><td>10</td><td>&nbsp;-</td><td>&nbsp;-</td><td>&nbsp;-</td><td></td><td>[3]</td></tr><tr><td>⋮</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr></tbody></table></figure>
<!-- /wp:table -->

<!-- wp:paragraph -->
<p>1.&nbsp;Qi-Jun Hong and Axel van de Walle,&nbsp;Solid-liquid coexistence in small systems: A statistical method to calculate melting temperatures.&nbsp;<em>Journal of Chemical Physics</em> <strong>139</strong> (9), 094114 (2013). [<a href="http://dx.doi.org/10.1063/1.4819792">DOI</a>]<br>2.&nbsp;Qi-Jun Hong, Sergey V. Ushakov, Alexandra Navrotsky, Axel van de Walle,&nbsp;Combined computational and experimental investigation of the refractory properties of La<sub>2</sub>Zr<sub>2</sub>O<sub>7</sub>.&nbsp;<em>Acta Materialia</em> <strong>84</strong>, 275-282 (2015). [<a href="http://dx.doi.org/10.1016/j.actamat.2014.10.026">DOI</a>]</p>
<!-- /wp:paragraph -->

<!-- wp:list {"ordered":true,"start":3} -->
<ol start="3" class="uds-list"><!-- wp:list-item -->
<li>Qi-Jun Hong and Axel van de Walle, Prediction of the material with highest melting temperature from quantum mechanics.&nbsp;<em>Physical Review B Rapid Communications</em>, <strong>92</strong>, 020104(R) (2015).&nbsp;[<a href="http://dx.doi.org/10.1103/PhysRevB.92.020104">DOI</a>]</li>
<!-- /wp:list-item -->

<!-- wp:list-item -->
<li>Ljubomir Miljacic, Steven Demers, Qi-Jun Hong and Axel van de Walle, Equation of state of solid, liquid and gaseous tantalum from first principles.&nbsp;<em>Calphad: Computer Coupling of Phase Diagrams and Thermochemistry</em>,&nbsp;<strong>51</strong>, 133-143 (2015). [<a href="http://dx.doi.org/10.1016/j.calphad.2015.08.005">DOI (open access)</a>].</li>
<!-- /wp:list-item --></ol>
<!-- /wp:list -->

<!-- wp:separator {"opacity":"css","className":"mt-12 is-style-copy-divider"} -->
<hr class="wp-block-separator has-css-opacity mt-12 is-style-copy-divider"/>
<!-- /wp:separator -->

<!-- /wp:wp-bootstrap-blocks/column -->
<!-- /wp:wp-bootstrap-blocks/row -->
<!-- /wp:wp-bootstrap-blocks/container -->
