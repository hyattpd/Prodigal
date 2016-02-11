/*******************************************************************************
    PRODIGAL (PROkaryotic DynamIc Programming Genefinding ALgorithm)
    Copyright (C) 2007-2016 University of Tennessee / UT-Battelle

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
*******************************************************************************/

#include "metagenomic.h"

/*******************************************************************************
  Initialize the metagenomic bins with the precalculated training files
  from the model organisms that best represent all of microbial Genbank.  
*******************************************************************************/
void initialize_metagenomic_bins(struct _metagenomic_bin *meta) {
  initialize_metagenome_0(meta[0].tinf);
  initialize_metagenome_1(meta[1].tinf);
  initialize_metagenome_2(meta[2].tinf);
  initialize_metagenome_3(meta[3].tinf);
  initialize_metagenome_4(meta[4].tinf);
  initialize_metagenome_5(meta[5].tinf);
  initialize_metagenome_6(meta[6].tinf);
  initialize_metagenome_7(meta[7].tinf);
  initialize_metagenome_8(meta[8].tinf);
  initialize_metagenome_9(meta[9].tinf);
  initialize_metagenome_10(meta[10].tinf);
  initialize_metagenome_11(meta[11].tinf);
  initialize_metagenome_12(meta[12].tinf);
  initialize_metagenome_13(meta[13].tinf);
  initialize_metagenome_14(meta[14].tinf);
  initialize_metagenome_15(meta[15].tinf);
  initialize_metagenome_16(meta[16].tinf);
  initialize_metagenome_17(meta[17].tinf);
  initialize_metagenome_18(meta[18].tinf);
  initialize_metagenome_19(meta[19].tinf);
  initialize_metagenome_20(meta[20].tinf);
  initialize_metagenome_21(meta[21].tinf);
  initialize_metagenome_22(meta[22].tinf);
  initialize_metagenome_23(meta[23].tinf);
  initialize_metagenome_24(meta[24].tinf);
  initialize_metagenome_25(meta[25].tinf);
  initialize_metagenome_26(meta[26].tinf);
  initialize_metagenome_27(meta[27].tinf);
  initialize_metagenome_28(meta[28].tinf);
  initialize_metagenome_29(meta[29].tinf);
  initialize_metagenome_30(meta[30].tinf);
  initialize_metagenome_31(meta[31].tinf);
  initialize_metagenome_32(meta[32].tinf);
  initialize_metagenome_33(meta[33].tinf);
  initialize_metagenome_34(meta[34].tinf);
  initialize_metagenome_35(meta[35].tinf);
  initialize_metagenome_36(meta[36].tinf);
  initialize_metagenome_37(meta[37].tinf);
  initialize_metagenome_38(meta[38].tinf);
  initialize_metagenome_39(meta[39].tinf);
  initialize_metagenome_40(meta[40].tinf);
  initialize_metagenome_41(meta[41].tinf);
  initialize_metagenome_42(meta[42].tinf);
  initialize_metagenome_43(meta[43].tinf);
  initialize_metagenome_44(meta[44].tinf);
  initialize_metagenome_45(meta[45].tinf);
  initialize_metagenome_46(meta[46].tinf);
  initialize_metagenome_47(meta[47].tinf);
  initialize_metagenome_48(meta[48].tinf);
  initialize_metagenome_49(meta[49].tinf);
  sprintf(meta[0].desc, "%d|%s|%s|%.1f|%d|%d", 0,
          "Mycoplasma_bovis_PG45",
          "B", 29.31, meta[0].tinf->trans_table, meta[0].tinf->uses_sd);
  sprintf(meta[1].desc, "%d|%s|%s|%.1f|%d|%d", 1,
          "Mycoplasma_pneumoniae_M129",
          "B", 40.01, meta[1].tinf->trans_table, meta[1].tinf->uses_sd);
  sprintf(meta[2].desc, "%d|%s|%s|%.1f|%d|%d", 2,
          "Mycoplasma_suis_Illinois",
          "B", 31.08, meta[2].tinf->trans_table, meta[2].tinf->uses_sd);
  sprintf(meta[3].desc, "%d|%s|%s|%.1f|%d|%d", 3,
          "Aeropyrum_pernix_K1",
          "A", 56.31, meta[3].tinf->trans_table, meta[3].tinf->uses_sd);
  sprintf(meta[4].desc, "%d|%s|%s|%.1f|%d|%d", 4,
          "Akkermansia_muciniphila_ATCC_BAA_835",
          "B", 55.76, meta[4].tinf->trans_table, meta[4].tinf->uses_sd);
  sprintf(meta[5].desc, "%d|%s|%s|%.1f|%d|%d", 5,
          "Anaplasma_marginale_Maries",
          "B", 49.76, meta[5].tinf->trans_table, meta[5].tinf->uses_sd);
  sprintf(meta[6].desc, "%d|%s|%s|%.1f|%d|%d", 6,
          "Anaplasma_phagocytophilum_HZ",
          "B", 41.64, meta[6].tinf->trans_table, meta[6].tinf->uses_sd);
  sprintf(meta[7].desc, "%d|%s|%s|%.1f|%d|%d", 7,
          "Archaeoglobus_fulgidus_DSM_4304",
          "A", 48.58, meta[7].tinf->trans_table, meta[7].tinf->uses_sd);
  sprintf(meta[8].desc, "%d|%s|%s|%.1f|%d|%d", 8,
          "Bacteroides_fragilis_NCTC_9343",
          "B", 43.19, meta[8].tinf->trans_table, meta[8].tinf->uses_sd);
  sprintf(meta[9].desc, "%d|%s|%s|%.1f|%d|%d", 9,
          "Brucella_canis_ATCC_23365",
          "B", 57.21, meta[9].tinf->trans_table, meta[9].tinf->uses_sd);
  sprintf(meta[10].desc, "%d|%s|%s|%.1f|%d|%d", 10,
          "Burkholderia_rhizoxinica_HKI_454",
          "B", 59.70, meta[10].tinf->trans_table, meta[10].tinf->uses_sd);
  sprintf(meta[11].desc, "%d|%s|%s|%.1f|%d|%d", 11,
          "Candidatus_Amoebophilus_asiaticus_5a2",
          "B", 35.05, meta[11].tinf->trans_table, meta[11].tinf->uses_sd);
  sprintf(meta[12].desc, "%d|%s|%s|%.1f|%d|%d", 12,
          "Candidatus_Korarchaeum_cryptofilum_OPF8",
          "A", 49.00, meta[12].tinf->trans_table, meta[12].tinf->uses_sd);
  sprintf(meta[13].desc, "%d|%s|%s|%.1f|%d|%d", 13,
          "Catenulispora_acidiphila_DSM_44928",
          "B", 69.77, meta[13].tinf->trans_table, meta[13].tinf->uses_sd);
  sprintf(meta[14].desc, "%d|%s|%s|%.1f|%d|%d", 14,
          "Cenarchaeum_symbiosum_B",
          "A", 57.19, meta[14].tinf->trans_table, meta[14].tinf->uses_sd);
  sprintf(meta[15].desc, "%d|%s|%s|%.1f|%d|%d", 15,
          "Chlorobium_phaeobacteroides_BS1",
          "B", 48.93, meta[15].tinf->trans_table, meta[15].tinf->uses_sd);
  sprintf(meta[16].desc, "%d|%s|%s|%.1f|%d|%d", 16,
          "Chlorobium_tepidum_TLS",
          "B", 56.53, meta[16].tinf->trans_table, meta[16].tinf->uses_sd);
  sprintf(meta[17].desc, "%d|%s|%s|%.1f|%d|%d", 17,
          "Desulfotomaculum_acetoxidans_DSM_771",
          "B", 41.55, meta[17].tinf->trans_table, meta[17].tinf->uses_sd);
  sprintf(meta[18].desc, "%d|%s|%s|%.1f|%d|%d", 18,
          "Desulfurococcus_kamchatkensis_1221n",
          "B", 45.34, meta[18].tinf->trans_table, meta[18].tinf->uses_sd);
  sprintf(meta[19].desc, "%d|%s|%s|%.1f|%d|%d", 19,
          "Erythrobacter_litoralis_HTCC2594",
          "B", 63.07, meta[19].tinf->trans_table, meta[19].tinf->uses_sd);
  sprintf(meta[20].desc, "%d|%s|%s|%.1f|%d|%d", 20,
          "Escherichia_coli_UMN026",
          "B", 50.72, meta[20].tinf->trans_table, meta[20].tinf->uses_sd);
  sprintf(meta[21].desc, "%d|%s|%s|%.1f|%d|%d", 21,
          "Haloquadratum_walsbyi_DSM_16790",
          "A", 47.86, meta[21].tinf->trans_table, meta[21].tinf->uses_sd);
  sprintf(meta[22].desc, "%d|%s|%s|%.1f|%d|%d", 22,
          "Halorubrum_lacusprofundi_ATCC_49239",
          "A", 57.14, meta[22].tinf->trans_table, meta[22].tinf->uses_sd);
  sprintf(meta[23].desc, "%d|%s|%s|%.1f|%d|%d", 23,
          "Hyperthermus_butylicus_DSM_5456",
          "A", 53.74, meta[23].tinf->trans_table, meta[23].tinf->uses_sd);
  sprintf(meta[24].desc, "%d|%s|%s|%.1f|%d|%d", 24,
          "Ignisphaera_aggregans_DSM_17230",
          "A", 35.69, meta[24].tinf->trans_table, meta[24].tinf->uses_sd);
  sprintf(meta[25].desc, "%d|%s|%s|%.1f|%d|%d", 25,
          "Marinobacter_aquaeolei_VT8",
          "B", 57.27, meta[25].tinf->trans_table, meta[25].tinf->uses_sd);
  sprintf(meta[26].desc, "%d|%s|%s|%.1f|%d|%d", 26,
          "Methanopyrus_kandleri_AV19",
          "A", 61.16, meta[26].tinf->trans_table, meta[26].tinf->uses_sd);
  sprintf(meta[27].desc, "%d|%s|%s|%.1f|%d|%d", 27,
          "Methanosphaerula_palustris_E1_9c",
          "A", 55.35, meta[27].tinf->trans_table, meta[27].tinf->uses_sd);
  sprintf(meta[28].desc, "%d|%s|%s|%.1f|%d|%d", 28,
          "Methanothermobacter_thermautotrophicus_Delta_H",
          "B", 49.54, meta[28].tinf->trans_table, meta[28].tinf->uses_sd);
  sprintf(meta[29].desc, "%d|%s|%s|%.1f|%d|%d", 29,
          "Methylacidiphilum_infernorum_V4",
          "B", 45.48, meta[29].tinf->trans_table, meta[29].tinf->uses_sd);
  sprintf(meta[30].desc, "%d|%s|%s|%.1f|%d|%d", 30,
          "Mycobacterium_leprae_TN",
          "B", 57.80, meta[30].tinf->trans_table, meta[30].tinf->uses_sd);
  sprintf(meta[31].desc, "%d|%s|%s|%.1f|%d|%d", 31,
          "Natrialba_magadii_ATCC_43099",
          "A", 61.42, meta[31].tinf->trans_table, meta[31].tinf->uses_sd);
  sprintf(meta[32].desc, "%d|%s|%s|%.1f|%d|%d", 32,
          "Orientia_tsutsugamushi_Boryong",
          "B", 30.53, meta[32].tinf->trans_table, meta[32].tinf->uses_sd);
  sprintf(meta[33].desc, "%d|%s|%s|%.1f|%d|%d", 33,
          "Pelotomaculum_thermopropionicum_SI",
          "B", 52.96, meta[33].tinf->trans_table, meta[33].tinf->uses_sd);
  sprintf(meta[34].desc, "%d|%s|%s|%.1f|%d|%d", 34,
          "Prochlorococcus_marinus_MIT_9313",
          "B", 50.74, meta[34].tinf->trans_table, meta[34].tinf->uses_sd);
  sprintf(meta[35].desc, "%d|%s|%s|%.1f|%d|%d", 35,
          "Pyrobaculum_aerophilum_IM2",
          "A", 51.36, meta[35].tinf->trans_table, meta[35].tinf->uses_sd);
  sprintf(meta[36].desc, "%d|%s|%s|%.1f|%d|%d", 36,
          "Ralstonia_solanacearum_PSI07",
          "B", 66.13, meta[36].tinf->trans_table, meta[36].tinf->uses_sd);
  sprintf(meta[37].desc, "%d|%s|%s|%.1f|%d|%d", 37,
          "Rhizobium_NGR234",
          "B", 58.49, meta[37].tinf->trans_table, meta[37].tinf->uses_sd);
  sprintf(meta[38].desc, "%d|%s|%s|%.1f|%d|%d", 38,
          "Rhodococcus_jostii_RHA1",
          "B", 65.05, meta[38].tinf->trans_table, meta[38].tinf->uses_sd);
  sprintf(meta[39].desc, "%d|%s|%s|%.1f|%d|%d", 39,
          "Rickettsia_conorii_Malish_7",
          "B", 32.44, meta[39].tinf->trans_table, meta[39].tinf->uses_sd);
  sprintf(meta[40].desc, "%d|%s|%s|%.1f|%d|%d", 40,
          "Rothia_dentocariosa_ATCC_17931",
          "B", 53.69, meta[40].tinf->trans_table, meta[40].tinf->uses_sd);
  sprintf(meta[41].desc, "%d|%s|%s|%.1f|%d|%d", 41,
          "Shigella_dysenteriae_Sd197",
          "B", 51.25, meta[41].tinf->trans_table, meta[41].tinf->uses_sd);
  sprintf(meta[42].desc, "%d|%s|%s|%.1f|%d|%d", 42,
          "Synechococcus_CC9605",
          "B", 59.22, meta[42].tinf->trans_table, meta[42].tinf->uses_sd);
  sprintf(meta[43].desc, "%d|%s|%s|%.1f|%d|%d", 43,
          "Synechococcus_JA_2_3B_a_2_13_",
          "B", 58.45, meta[43].tinf->trans_table, meta[43].tinf->uses_sd);
  sprintf(meta[44].desc, "%d|%s|%s|%.1f|%d|%d", 44,
          "Thermoplasma_volcanium_GSS1",
          "A", 39.92, meta[44].tinf->trans_table, meta[44].tinf->uses_sd);
  sprintf(meta[45].desc, "%d|%s|%s|%.1f|%d|%d", 45,
          "Treponema_pallidum_Nichols",
          "B", 52.77, meta[45].tinf->trans_table, meta[45].tinf->uses_sd);
  sprintf(meta[46].desc, "%d|%s|%s|%.1f|%d|%d", 46,
          "Tropheryma_whipplei_TW08_27",
          "B", 46.31, meta[46].tinf->trans_table, meta[46].tinf->uses_sd);
  sprintf(meta[47].desc, "%d|%s|%s|%.1f|%d|%d", 47,
          "Xenorhabdus_nematophila_ATCC_19061",
          "B", 44.15, meta[47].tinf->trans_table, meta[47].tinf->uses_sd);
  sprintf(meta[48].desc, "%d|%s|%s|%.1f|%d|%d", 48,
          "Xylella_fastidiosa_Temecula1",
          "B", 51.78, meta[48].tinf->trans_table, meta[48].tinf->uses_sd);
  sprintf(meta[49].desc, "%d|%s|%s|%.1f|%d|%d", 49,
          "_Nostoc_azollae__0708",
          "B", 38.45, meta[49].tinf->trans_table, meta[49].tinf->uses_sd);
}
