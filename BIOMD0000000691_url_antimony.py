import tellurium as te

r = te.loada('''
// Created by libAntimony v2.12.0.3
function Constant_flux__irreversible(v)
  v;
end

Constant_flux__irreversible is "Constant flux (irreversible)"

function function_for_v2(k2, S2, N1)
  k2*S2*N1;
end

function_for_v2 is "function for v2"

function function_for_d_dt_S4_ex(phi, J1, J2)
  (phi/2)*(J1 + J2);
end

function_for_d_dt_S4_ex is "function for d/dt(S4_ex)"

function function_for_v1(k1, S1, A3, K_I, q)
  k1*S1*A3*(1 + (A3/K_I)^q)^-1;
end

function_for_v1 is "function for v1"

function function_for_v3(k3, S3, A2)
  k3*S3*A2;
end

function_for_v3 is "function for v3"


model *Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast()

  // Compartments and Species:
  compartment Cell_1, Cell_2, Compartment_;
  species S1__Cell_1_ in Cell_1, S1__Cell_2_ in Cell_2, S2__Cell_1_ in Cell_1;
  species S2__Cell_2_ in Cell_2, S3__Cell_1_ in Cell_1, S3__Cell_2_ in Cell_2;
  species S4__Cell_1_ in Cell_1, S4__Cell_2_ in Cell_2, N2__Cell_1_ in Cell_1;
  species N2__Cell_2_ in Cell_2, A3__Cell_1_ in Cell_1, A3__Cell_2_ in Cell_2;
  species S4_ex in Compartment_, $A in Compartment_, $N in Compartment_, $N1__Cell_1_ in Cell_1;
  species $N1__Cell_2_ in Cell_2, $A2__Cell_1_ in Cell_1, $A2__Cell_2_ in Cell_2;

  // Assignment Rules:
  N1__Cell_1_ := N - N2__Cell_1_;
  N1__Cell_2_ := N - N2__Cell_2_;
  A2__Cell_1_ := A - A3__Cell_1_;
  A2__Cell_2_ := A - A3__Cell_2_;
  J_cell_1 := kappa*(S4__Cell_1_ - S4_ex);
  j_cell_2 := kappa*(S4__Cell_2_ - S4_ex);

  // Reactions:
  v1_cell_1: S1__Cell_1_ + 2 A3__Cell_1_ => 2 S2__Cell_1_; Cell_1*function_for_v1(k1, S1__Cell_1_, A3__Cell_1_, K_I, q);
  v1_cell_2: S1__Cell_2_ + 2 A3__Cell_2_ => 2 S2__Cell_2_; Cell_2*function_for_v1(k1, S1__Cell_2_, A3__Cell_2_, K_I, q);
  v2_cell_1: S2__Cell_1_ + $N1__Cell_1_ => S3__Cell_1_ + N2__Cell_1_; Cell_1*function_for_v2(k2, S2__Cell_1_, N1__Cell_1_);
  v2_cell_2: S2__Cell_2_ + $N1__Cell_2_ => S3__Cell_2_ + N2__Cell_2_; Cell_2*function_for_v2(k2, S2__Cell_2_, N1__Cell_2_);
  v3_cell_1: S3__Cell_1_ + $A2__Cell_1_ => S4__Cell_1_ + 2 A3__Cell_1_; Cell_1*function_for_v3(k3, S3__Cell_1_, A2__Cell_1_);
  v3_cell_2: S3__Cell_2_ + $A2__Cell_2_ => S4__Cell_2_ + 2 A3__Cell_2_; Cell_2*function_for_v3(k3, S3__Cell_2_, A2__Cell_2_);
  v4_cell_1: S4__Cell_1_ + N2__Cell_1_ => ; Cell_1*k4*S4__Cell_1_*N2__Cell_1_;
  v4_cell_2: S4__Cell_2_ + N2__Cell_2_ => ; Cell_2*k4*S4__Cell_2_*N2__Cell_2_;
  v5_cell_1: A3__Cell_1_ => ; Cell_1*k5*A3__Cell_1_;
  v5_cell_2: A3__Cell_2_ => ; Cell_2*k5*A3__Cell_2_;
  v6_cell_1: S2__Cell_1_ + N2__Cell_1_ => ; Cell_1*k6*S2__Cell_1_*N2__Cell_1_;
  v6_cell_2: S2__Cell_2_ + N2__Cell_2_ => ; Cell_2*k6*S2__Cell_2_*N2__Cell_2_;
  v7: S4_ex => ; Compartment_*k*S4_ex;
  S1_cell_1_glucose_influx:  => S1__Cell_1_; Cell_1*Constant_flux__irreversible(J0);
  S1_cell_2_glucose_influx:  => S1__Cell_2_; Cell_2*Constant_flux__irreversible(J0);
  S4_cell_1_export: S4__Cell_1_ => ; Cell_1*Constant_flux__irreversible(J_cell_1);
  S4_cell_2_export: S4__Cell_2_ => ; Cell_2*Constant_flux__irreversible(j_cell_2);
  S4_ex_import:  => S4_ex; Compartment_*function_for_d_dt_S4_ex(phi, J_cell_1, j_cell_2);

  // Species initializations:
  S1__Cell_1_ = 5.8;
  S1__Cell_2_ = 2.9;
  S2__Cell_1_ = 0.9;
  S2__Cell_2_ = 0.45;
  S3__Cell_1_ = 0.2;
  S3__Cell_2_ = 0.1;
  S4__Cell_1_ = 0.2;
  S4__Cell_2_ = 0.1;
  N2__Cell_1_ = 0.1;
  N2__Cell_2_ = 0.05;
  A3__Cell_1_ = 3.2;
  A3__Cell_2_ = 0.2;
  S4_ex = 0.1;
  A = 4;
  N = 1;

  // Compartment initializations:
  Cell_1 = 1;
  Cell_2 = 1;
  Compartment_ = 1;

  // Variable initializations:
  k1 = 100;
  K_I = 0.52;
  q = 4;
  k2 = 6;
  k3 = 16;
  k4 = 100;
  k5 = 1.28;
  k6 = 12;
  k = 1.5;
  J0 = 3;
  kappa = 13;
  phi = 0.1;

  // Other declarations:
  var J_cell_1, j_cell_2;
  const Cell_1, Cell_2, Compartment_, k1, K_I, q, k2, k3, k4, k5, k6, k, J0;
  const kappa, phi;

  // Unit definitions:
  unit volume = 1e-3 litre;
  unit substance = 1e-3 mole;

  // Display Names:
  Cell_1 is "Cell 1";
  Cell_2 is "Cell 2";
  S1__Cell_1_ is "S1";
  S1__Cell_2_ is "S1";
  S2__Cell_1_ is "S2";
  S2__Cell_2_ is "S2";
  S3__Cell_1_ is "S3";
  S3__Cell_2_ is "S3";
  S4__Cell_1_ is "S4";
  S4__Cell_2_ is "S4";
  N2__Cell_1_ is "N2";
  N2__Cell_2_ is "N2";
  A3__Cell_1_ is "A3";
  A3__Cell_2_ is "A3";
  N1__Cell_1_ is "N1";
  N1__Cell_2_ is "N1";
  A2__Cell_1_ is "A2";
  A2__Cell_2_ is "A2";
  S1_cell_1_glucose_influx is "S1_cell_1 glucose influx";
  S1_cell_2_glucose_influx is "S1_cell_2 glucose influx";
  S4_cell_1_export is "S4_cell_1 export";
  S4_cell_2_export is "S4_cell_2 export";
  S4_ex_import is "S4_ex import";

  // CV terms:
  Cell_1 hypernym "http://identifiers.org/taxonomy/4930"
  Cell_2 hypernym "http://identifiers.org/taxonomy/4930"
  S1__Cell_1_ hypernym "http://identifiers.org/chebi/CHEBI:17234"
  S1__Cell_1_ hypernym "http://identifiers.org/kegg.compound/C00031"
  S1__Cell_2_ hypernym "http://identifiers.org/chebi/CHEBI:17234"
  S1__Cell_2_ hypernym "http://identifiers.org/kegg.compound/C00031"
  S2__Cell_1_ part "http://identifiers.org/kegg.compound/C00661"
  S2__Cell_1_ part "http://identifiers.org/kegg.compound/C00111"
  S2__Cell_2_ part "http://identifiers.org/kegg.compound/C00661"
  S2__Cell_2_ part "http://identifiers.org/kegg.compound/C00111"
  S3__Cell_1_ hypernym "http://identifiers.org/chebi/CHEBI:89363"
  S3__Cell_1_ hypernym "http://identifiers.org/kegg.compound/C00236"
  S3__Cell_2_ hypernym "http://identifiers.org/chebi/CHEBI:89363"
  S3__Cell_2_ hypernym "http://identifiers.org/kegg.compound/C00236"
  S4__Cell_1_ part "http://identifiers.org/kegg.compound/C00022"
  S4__Cell_1_ part "http://identifiers.org/kegg.compound/C00084"
  S4__Cell_2_ part "http://identifiers.org/kegg.compound/C00084"
  S4__Cell_2_ part "http://identifiers.org/kegg.compound/C00022"
  N2__Cell_1_ hypernym "http://identifiers.org/kegg.compound/C00004"
  N2__Cell_1_ hypernym "http://identifiers.org/chebi/CHEBI:16908"
  N2__Cell_2_ hypernym "http://identifiers.org/chebi/CHEBI:16908"
  N2__Cell_2_ hypernym "http://identifiers.org/kegg.compound/C00004"
  A3__Cell_1_ hypernym "http://identifiers.org/chebi/CHEBI:15422"
  A3__Cell_1_ hypernym "http://identifiers.org/kegg.compound/C00002"
  A3__Cell_2_ hypernym "http://identifiers.org/chebi/CHEBI:15422"
  A3__Cell_2_ hypernym "http://identifiers.org/kegg.compound/C00002"
  A part "http://identifiers.org/kegg.compound/C00002"
  A part "http://identifiers.org/kegg.compound/C00008"
  N part "http://identifiers.org/kegg.compound/C00004"
  N part "http://identifiers.org/kegg.compound/C00003"
  N1__Cell_1_ hypernym "http://identifiers.org/chebi/CHEBI:15846"
  N1__Cell_1_ hypernym "http://identifiers.org/kegg.compound/C00003"
  N1__Cell_2_ hypernym "http://identifiers.org/kegg.compound/C00003"
  N1__Cell_2_ hypernym "http://identifiers.org/chebi/CHEBI:15846"
  A2__Cell_1_ hypernym "http://identifiers.org/chebi/CHEBI:16761"
  A2__Cell_1_ hypernym "http://identifiers.org/kegg.compound/C00008"
  A2__Cell_2_ hypernym "http://identifiers.org/chebi/CHEBI:16761"
  A2__Cell_2_ hypernym "http://identifiers.org/kegg.compound/C00008"
  K_I encodement "http://identifiers.org/sbo/SBO:0000169"
  kappa encodement "http://identifiers.org/go/GO:0090559"
  v1_cell_1 part "http://identifiers.org/kegg.reaction/R01015"
  v1_cell_1 part "http://identifiers.org/kegg.reaction/R01070"
  v1_cell_1 part "http://identifiers.org/kegg.reaction/R02740"
  v1_cell_1 part "http://identifiers.org/kegg.reaction/R04779"
  v1_cell_1 part "http://identifiers.org/kegg.reaction/R01786"
  v1_cell_2 part "http://identifiers.org/kegg.reaction/R01070"
  v1_cell_2 part "http://identifiers.org/kegg.reaction/R01015"
  v1_cell_2 part "http://identifiers.org/kegg.reaction/R02740"
  v1_cell_2 part "http://identifiers.org/kegg.reaction/R04779"
  v1_cell_2 part "http://identifiers.org/kegg.reaction/R01786"
  v2_cell_1 hypernym "http://identifiers.org/kegg.reaction/R01061"
  v2_cell_2 hypernym "http://identifiers.org/kegg.reaction/R01061"
  v3_cell_1 part "http://identifiers.org/kegg.reaction/R00200"
  v3_cell_1 part "http://identifiers.org/kegg.reaction/R00658"
  v3_cell_1 part "http://identifiers.org/kegg.reaction/R01512"
  v3_cell_1 part "http://identifiers.org/kegg.reaction/R01518"
  v3_cell_2 part "http://identifiers.org/kegg.reaction/R01512"
  v3_cell_2 part "http://identifiers.org/kegg.reaction/R00658"
  v3_cell_2 part "http://identifiers.org/kegg.reaction/R00200"
  v3_cell_2 part "http://identifiers.org/kegg.reaction/R01518"
  v4_cell_1 hypernym "http://identifiers.org/kegg.reaction/R00754"
  v4_cell_2 hypernym "http://identifiers.org/kegg.reaction/R00754"
  v6_cell_1 hypernym "http://identifiers.org/kegg.reaction/R01036"
  v6_cell_2 hypernym "http://identifiers.org/kegg.reaction/R01036"
  v7 hypernym "http://identifiers.org/sbo/SBO:0000179"
  S1_cell_1_glucose_influx hypernym "http://identifiers.org/go/GO:1904659"
  S1_cell_2_glucose_influx hypernym "http://identifiers.org/go/GO:1904659"
  S4_cell_1_export hypernym "http://identifiers.org/go/GO:1901475"
  S4_cell_2_export hypernym "http://identifiers.org/go/GO:1901475"
  S4_ex_import hypernym "http://identifiers.org/go/GO:1901475"
end

Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast is "Wolf2000 - Cellular interaction on glycolytic oscillations in yeast"

Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast model_entity_is "http://identifiers.org/biomodels.db/MODEL1006230022"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast description "http://identifiers.org/pubmed/10702114"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast model_entity_is "http://identifiers.org/biomodels.db/BIOMD0000000691"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast hypernym "http://identifiers.org/reactome/R-HSA-70171"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast hypernym "http://identifiers.org/go/GO:0006096"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast hypernym "http://identifiers.org/kegg.pathway/map00010"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast property "http://identifiers.org/go/GO:0009061"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast taxon "http://identifiers.org/taxonomy/4930"
Wolf2000___Cellular_interaction_on_glycolytic_oscillations_in_yeast description "http://identifiers.org/pubmed/10702114"
''')

m = r.simulate (0, 10, 100)
r.plot()