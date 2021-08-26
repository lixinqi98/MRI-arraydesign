## Project Record

Record the project based on the timeline.

1. Re-draw the UNIC 51, UNIC 68 and UNIC 99 based on the previous design. Detail figures can be seen in the excel.

   + UNIC 68: original code `BiotSavart_MultiCoil_JH_2_5_2_1.m`.
   + UNIC 53: based on `BiotSavart_MultiCoil_JH_2_5_2_1.m`, remove the loops marked as black.
   + UNIC 99: `MonaBiotSavart_MultiCoil_JH_2_5_2_1.m`
   + more coils(add small loops inside the big loops, not used): `Mona_v2_BiotSavart_MultiCoil_JH_2_5_2_1.m`

2. Create the **SH7** based on the formula and get the 7T original field maps from @Chris.

   + `SH_order_20210527.m`

3. Calculate the **3T** whole brain shimming results of SH2-7, iPRES 32, UNIC 52, UNIC 68, and UNIC 99. Compare shimming effects with the previous one and check the field map before&after.

   + use the same brain&coil position as @Carissa, `offsety = -0.01; offsetz = 0.02`
   + UNIC: `MainWholeBrainShimming_UNIC_3T.m`
   + SH:`MainWholeBrainShimming_SHorder_3T.m`

4. Apply the multi-edge cutting to the brains and compare them with results in 3. Use the field map with multi-edge cutting in all following calculation.

   + `B0_multi_Edgecut.m`

5. Calculate the **7T** whole brain shimming results of SH2-7, iPRES 32, UNIC 51, UNIC 68, and UNIC 99. 

   + **The most important difference between 3T and 7T is that in 3T we have to scale down 2pi and in 7T no need to do this.** (How we find this? After checking the detailed field map of each figure and each slice, found that the average std of each slice is strange and different from the common one.) 
   + The loops' current limit should be differnet in 3T and 7T. (add a detail table here)
   + **7T field maps have different original sizes and they have to be aligned to the same size**. We find the largest x,y,z of all brains and add `nan` padding to each brain. 
   + UNIC: `MainWholeBrainShimming_UNIC_7T.m`
     + **Use the correct brain position** The relative position of brain and coils has a great influence on the shimming results. We tested 4 different settings as shown in excel and finally used `offsety = -0.020; offsetz = 0.03;`
   + SH: `MainWholeBrainShimming_SHorder_7T.m `

6. After making sure the correction of 3T&7T whole brain shimming, we do the **slab shimming**. Each slab consists of 5 slices and each slab's shimming effect is calculated independently.( add a table of slab details)

   + UNIC: `MainSlabShimming_UNIC_3T.m`, `MainSlabShimming_UNIC_7T.m`
   + SH: `MainSlabShimming_SHorder_3T.m`, `MainSlabShimming_SHorder_7T.m`

7. Due to the redundant coil design in UNIC 52/68/99, we apply the coil reduction. The coil reduction takes the advantage of the idea of PCA. The reserved coils' contribution to important components in field map was sorted and the least important loops are moved in iterations. 

   + Integrated in the `Main*` codes, usually in `STEP 5`

8. In first version of coil reduction, coils' weights are calculated based on whole brain shimming, so we can not guarantee the effect of slab shimming. In the next version, the slab shimming effect was taken into consideration. Each loop's weight was calculated based on the ensemble of whole brain and slab shimming. Like here:

   ```
   coil weight = 0.5 * `wholegrain_weight` + 0.5 * (`slab1_weight`* 1/6 + `slab2_weight`* 1/6 + `slab3_weight`* 1/6 + `slab4_weight`* 1/6+ `slab5_weight`* 1/6 + `slab6_weight`* 1/6)
   ```

   , which can be modified in `ChooseCoil_WorldCordf_v2.m` 

9. The coil reduction is based on the dataset of 16 brains and to test the robust of our reseved coil design, we have to test the whole brain and slab shimming effects on the test set. 8 extra brains are randomly picked from the HCP dataset. In conclusion, for 3T&7T, we use **16 subjects train set and 8 subjects test set**.
