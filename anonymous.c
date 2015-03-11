/******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2015 University of Tennessee / UT-Battelle

    Code Author:  Doug Hyatt

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#include "anonymous.h"

/******************************************************************************
  Given GC content, return minimum and maximum GC boundaries within
  which to try anonymous training files.  The formula for these
  bounds was determined through linear regression examining a million
  gene fragments of various sizes.
******************************************************************************/
void get_gc_bounds(double *bound, double gc)
{
  bound[0] = 0.88495 * gc - 0.0102337;
  if (bound[0] > 0.65)
  {
    bound[0] = 0.65;
  }
  bound[1] = 0.86596 * gc + .1131991;
  if (bound[1] < 0.35)
  {
    bound[1] = 0.35;
  }
}

/******************************************************************************
  Initialize the preset genome bins with the precalculated training files
  from the model organisms that best represent all of microbial Genbank.
******************************************************************************/
void initialize_preset_genome_bins(struct _preset_genome_bin *presets)
{
  initialize_preset_genome_0(presets[0].data);
  initialize_preset_genome_1(presets[1].data);
  initialize_preset_genome_2(presets[2].data);
  initialize_preset_genome_3(presets[3].data);
  initialize_preset_genome_4(presets[4].data);
  initialize_preset_genome_5(presets[5].data);
  initialize_preset_genome_6(presets[6].data);
  initialize_preset_genome_7(presets[7].data);
  initialize_preset_genome_8(presets[8].data);
  initialize_preset_genome_9(presets[9].data);
  initialize_preset_genome_10(presets[10].data);
  initialize_preset_genome_11(presets[11].data);
  initialize_preset_genome_12(presets[12].data);
  initialize_preset_genome_13(presets[13].data);
  initialize_preset_genome_14(presets[14].data);
  initialize_preset_genome_15(presets[15].data);
  initialize_preset_genome_16(presets[16].data);
  initialize_preset_genome_17(presets[17].data);
  initialize_preset_genome_18(presets[18].data);
  initialize_preset_genome_19(presets[19].data);
  initialize_preset_genome_20(presets[20].data);
  initialize_preset_genome_21(presets[21].data);
  initialize_preset_genome_22(presets[22].data);
  initialize_preset_genome_23(presets[23].data);
  initialize_preset_genome_24(presets[24].data);
  initialize_preset_genome_25(presets[25].data);
  initialize_preset_genome_26(presets[26].data);
  initialize_preset_genome_27(presets[27].data);
  initialize_preset_genome_28(presets[28].data);
  initialize_preset_genome_29(presets[29].data);
  initialize_preset_genome_30(presets[30].data);
  initialize_preset_genome_31(presets[31].data);
  initialize_preset_genome_32(presets[32].data);
  initialize_preset_genome_33(presets[33].data);
  initialize_preset_genome_34(presets[34].data);
  initialize_preset_genome_35(presets[35].data);
  initialize_preset_genome_36(presets[36].data);
  initialize_preset_genome_37(presets[37].data);
  initialize_preset_genome_38(presets[38].data);
  initialize_preset_genome_39(presets[39].data);
  initialize_preset_genome_40(presets[40].data);
  initialize_preset_genome_41(presets[41].data);
  initialize_preset_genome_42(presets[42].data);
  initialize_preset_genome_43(presets[43].data);
  initialize_preset_genome_44(presets[44].data);
  initialize_preset_genome_45(presets[45].data);
  initialize_preset_genome_46(presets[46].data);
  initialize_preset_genome_47(presets[47].data);
  initialize_preset_genome_48(presets[48].data);
  initialize_preset_genome_49(presets[49].data);
  sprintf(presets[0].desc, "%d|%s|%s|%.1f|%d|%d", 0,
          "Mycoplasma_bovis_PG45",
          "B", 29.31, presets[0].data->trans_table, presets[0].data->uses_sd);
  sprintf(presets[1].desc, "%d|%s|%s|%.1f|%d|%d", 1,
          "Mycoplasma_pneumoniae_M129",
          "B", 40.01, presets[1].data->trans_table, presets[1].data->uses_sd);
  sprintf(presets[2].desc, "%d|%s|%s|%.1f|%d|%d", 2,
          "Mycoplasma_suis_Illinois",
          "B", 31.08, presets[2].data->trans_table, presets[2].data->uses_sd);
  sprintf(presets[3].desc, "%d|%s|%s|%.1f|%d|%d", 3,
          "Aeropyrum_pernix_K1",
          "A", 56.31, presets[3].data->trans_table, presets[3].data->uses_sd);
  sprintf(presets[4].desc, "%d|%s|%s|%.1f|%d|%d", 4,
          "Akkermansia_muciniphila_ATCC_BAA_835",
          "B", 55.76, presets[4].data->trans_table, presets[4].data->uses_sd);
  sprintf(presets[5].desc, "%d|%s|%s|%.1f|%d|%d", 5,
          "Anaplasma_marginale_Maries",
          "B", 49.76, presets[5].data->trans_table, presets[5].data->uses_sd);
  sprintf(presets[6].desc, "%d|%s|%s|%.1f|%d|%d", 6,
          "Anaplasma_phagocytophilum_HZ",
          "B", 41.64, presets[6].data->trans_table, presets[6].data->uses_sd);
  sprintf(presets[7].desc, "%d|%s|%s|%.1f|%d|%d", 7,
          "Archaeoglobus_fulgidus_DSM_4304",
          "A", 48.58, presets[7].data->trans_table, presets[7].data->uses_sd);
  sprintf(presets[8].desc, "%d|%s|%s|%.1f|%d|%d", 8,
          "Bacteroides_fragilis_NCTC_9343",
          "B", 43.19, presets[8].data->trans_table, presets[8].data->uses_sd);
  sprintf(presets[9].desc, "%d|%s|%s|%.1f|%d|%d", 9,
          "Brucella_canis_ATCC_23365",
          "B", 57.21, presets[9].data->trans_table, presets[9].data->uses_sd);
  sprintf(presets[10].desc, "%d|%s|%s|%.1f|%d|%d", 10,
          "Burkholderia_rhizoxinica_HKI_454",
          "B",59.70, presets[10].data->trans_table, presets[10].data->uses_sd);
  sprintf(presets[11].desc, "%d|%s|%s|%.1f|%d|%d", 11,
          "Candidatus_Amoebophilus_asiaticus_5a2",
          "B",35.05, presets[11].data->trans_table, presets[11].data->uses_sd);
  sprintf(presets[12].desc, "%d|%s|%s|%.1f|%d|%d", 12,
          "Candidatus_Korarchaeum_cryptofilum_OPF8",
          "A",49.00, presets[12].data->trans_table, presets[12].data->uses_sd);
  sprintf(presets[13].desc, "%d|%s|%s|%.1f|%d|%d", 13,
          "Catenulispora_acidiphila_DSM_44928",
          "B",69.77, presets[13].data->trans_table, presets[13].data->uses_sd);
  sprintf(presets[14].desc, "%d|%s|%s|%.1f|%d|%d", 14,
          "Cenarchaeum_symbiosum_B",
          "A",57.19, presets[14].data->trans_table, presets[14].data->uses_sd);
  sprintf(presets[15].desc, "%d|%s|%s|%.1f|%d|%d", 15,
          "Chlorobium_phaeobacteroides_BS1",
          "B",48.93, presets[15].data->trans_table, presets[15].data->uses_sd);
  sprintf(presets[16].desc, "%d|%s|%s|%.1f|%d|%d", 16,
          "Chlorobium_tepidum_TLS",
          "B",56.53, presets[16].data->trans_table, presets[16].data->uses_sd);
  sprintf(presets[17].desc, "%d|%s|%s|%.1f|%d|%d", 17,
          "Desulfotomaculum_acetoxidans_DSM_771",
          "B",41.55, presets[17].data->trans_table, presets[17].data->uses_sd);
  sprintf(presets[18].desc, "%d|%s|%s|%.1f|%d|%d", 18,
          "Desulfurococcus_kamchatkensis_1221n",
          "B",45.34, presets[18].data->trans_table, presets[18].data->uses_sd);
  sprintf(presets[19].desc, "%d|%s|%s|%.1f|%d|%d", 19,
          "Erythrobacter_litoralis_HTCC2594",
          "B",63.07, presets[19].data->trans_table, presets[19].data->uses_sd);
  sprintf(presets[20].desc, "%d|%s|%s|%.1f|%d|%d", 20,
          "Escherichia_coli_UMN026",
          "B",50.72, presets[20].data->trans_table, presets[20].data->uses_sd);
  sprintf(presets[21].desc, "%d|%s|%s|%.1f|%d|%d", 21,
          "Haloquadratum_walsbyi_DSM_16790",
          "A",47.86, presets[21].data->trans_table, presets[21].data->uses_sd);
  sprintf(presets[22].desc, "%d|%s|%s|%.1f|%d|%d", 22,
          "Halorubrum_lacusprofundi_ATCC_49239",
          "A",57.14, presets[22].data->trans_table, presets[22].data->uses_sd);
  sprintf(presets[23].desc, "%d|%s|%s|%.1f|%d|%d", 23,
          "Hyperthermus_butylicus_DSM_5456",
          "A",53.74, presets[23].data->trans_table, presets[23].data->uses_sd);
  sprintf(presets[24].desc, "%d|%s|%s|%.1f|%d|%d", 24,
          "Ignisphaera_aggregans_DSM_17230",
          "A",35.69, presets[24].data->trans_table, presets[24].data->uses_sd);
  sprintf(presets[25].desc, "%d|%s|%s|%.1f|%d|%d", 25,
          "Marinobacter_aquaeolei_VT8",
          "B",57.27, presets[25].data->trans_table, presets[25].data->uses_sd);
  sprintf(presets[26].desc, "%d|%s|%s|%.1f|%d|%d", 26,
          "Methanopyrus_kandleri_AV19",
          "A",61.16, presets[26].data->trans_table, presets[26].data->uses_sd);
  sprintf(presets[27].desc, "%d|%s|%s|%.1f|%d|%d", 27,
          "Methanosphaerula_palustris_E1_9c",
          "A",55.35, presets[27].data->trans_table, presets[27].data->uses_sd);
  sprintf(presets[28].desc, "%d|%s|%s|%.1f|%d|%d", 28,
          "Methanothermobacter_thermautotrophicus_Delta_H",
          "B",49.54, presets[28].data->trans_table, presets[28].data->uses_sd);
  sprintf(presets[29].desc, "%d|%s|%s|%.1f|%d|%d", 29,
          "Methylacidiphilum_infernorum_V4",
          "B",45.48, presets[29].data->trans_table, presets[29].data->uses_sd);
  sprintf(presets[30].desc, "%d|%s|%s|%.1f|%d|%d", 30,
          "Mycobacterium_leprae_TN",
          "B",57.80, presets[30].data->trans_table, presets[30].data->uses_sd);
  sprintf(presets[31].desc, "%d|%s|%s|%.1f|%d|%d", 31,
          "Natrialba_magadii_ATCC_43099",
          "A",61.42, presets[31].data->trans_table, presets[31].data->uses_sd);
  sprintf(presets[32].desc, "%d|%s|%s|%.1f|%d|%d", 32,
          "Orientia_tsutsugamushi_Boryong",
          "B",30.53, presets[32].data->trans_table, presets[32].data->uses_sd);
  sprintf(presets[33].desc, "%d|%s|%s|%.1f|%d|%d", 33,
          "Pelotomaculum_thermopropionicum_SI",
          "B",52.96, presets[33].data->trans_table, presets[33].data->uses_sd);
  sprintf(presets[34].desc, "%d|%s|%s|%.1f|%d|%d", 34,
          "Prochlorococcus_marinus_MIT_9313",
          "B",50.74, presets[34].data->trans_table, presets[34].data->uses_sd);
  sprintf(presets[35].desc, "%d|%s|%s|%.1f|%d|%d", 35,
          "Pyrobaculum_aerophilum_IM2",
          "A",51.36, presets[35].data->trans_table, presets[35].data->uses_sd);
  sprintf(presets[36].desc, "%d|%s|%s|%.1f|%d|%d", 36,
          "Ralstonia_solanacearum_PSI07",
          "B",66.13, presets[36].data->trans_table, presets[36].data->uses_sd);
  sprintf(presets[37].desc, "%d|%s|%s|%.1f|%d|%d", 37,
          "Rhizobium_NGR234",
          "B",58.49, presets[37].data->trans_table, presets[37].data->uses_sd);
  sprintf(presets[38].desc, "%d|%s|%s|%.1f|%d|%d", 38,
          "Rhodococcus_jostii_RHA1",
          "B",65.05, presets[38].data->trans_table, presets[38].data->uses_sd);
  sprintf(presets[39].desc, "%d|%s|%s|%.1f|%d|%d", 39,
          "Rickettsia_conorii_Malish_7",
          "B",32.44, presets[39].data->trans_table, presets[39].data->uses_sd);
  sprintf(presets[40].desc, "%d|%s|%s|%.1f|%d|%d", 40,
          "Rothia_dentocariosa_ATCC_17931",
          "B",53.69, presets[40].data->trans_table, presets[40].data->uses_sd);
  sprintf(presets[41].desc, "%d|%s|%s|%.1f|%d|%d", 41,
          "Shigella_dysenteriae_Sd197",
          "B",51.25, presets[41].data->trans_table, presets[41].data->uses_sd);
  sprintf(presets[42].desc, "%d|%s|%s|%.1f|%d|%d", 42,
          "Synechococcus_CC9605",
          "B",59.22, presets[42].data->trans_table, presets[42].data->uses_sd);
  sprintf(presets[43].desc, "%d|%s|%s|%.1f|%d|%d", 43,
          "Synechococcus_JA_2_3B_a_2_13_",
          "B",58.45, presets[43].data->trans_table, presets[43].data->uses_sd);
  sprintf(presets[44].desc, "%d|%s|%s|%.1f|%d|%d", 44,
          "Thermoplasma_volcanium_GSS1",
          "A",39.92, presets[44].data->trans_table, presets[44].data->uses_sd);
  sprintf(presets[45].desc, "%d|%s|%s|%.1f|%d|%d", 45,
          "Treponema_pallidum_Nichols",
          "B",52.77, presets[45].data->trans_table, presets[45].data->uses_sd);
  sprintf(presets[46].desc, "%d|%s|%s|%.1f|%d|%d", 46,
          "Tropheryma_whipplei_TW08_27",
          "B",46.31, presets[46].data->trans_table, presets[46].data->uses_sd);
  sprintf(presets[47].desc, "%d|%s|%s|%.1f|%d|%d", 47,
          "Xenorhabdus_nematophila_ATCC_19061",
          "B",44.15, presets[47].data->trans_table, presets[47].data->uses_sd);
  sprintf(presets[48].desc, "%d|%s|%s|%.1f|%d|%d", 48,
          "Xylella_fastidiosa_Temecula1",
          "B",51.78, presets[48].data->trans_table, presets[48].data->uses_sd);
  sprintf(presets[49].desc, "%d|%s|%s|%.1f|%d|%d", 49,
          "_Nostoc_azollae__0708",
          "B",38.45, presets[49].data->trans_table, presets[49].data->uses_sd);
}
