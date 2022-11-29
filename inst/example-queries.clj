;;; Get me all the nanostring count measurements for the Roh2017 dataset

(d/q '[:find ?sample-id ?gene-hugo ?count
       :in $ ?dataset
       :where
       [?d :dataset/name ?dataset]
       [?d :dataset/assays ?a]
       [?a :assay/technology :technology/nanostring]
       [?a :assay/measurement-sets ?e]
       [?e :measurement-set/measurements ?m]
       [?m :measurement/nanostring-count ?count]
       [?m :measurement/gene ?g]
       [?g :gene/hugo ?gene-hugo]
       [?m :measurement/sample ?s]
       [?s :sample/id ?sample-id]]
     db "Roh2017")


;;; Get me all the variants for the Roh2017 dataset

(d/q '[:find ?sample-id ?gene-hugo ?type
       :in $ ?dataset
       :where
       [?d :dataset/name ?dataset]
       [?d :dataset/assays ?a]
       [?a :assay/technology :technology/WES]
       [?a :assay/measurement-sets ?e]
       [?e :measurement-set/measurements ?m]
       [?m :measurement/variant ?v]
       [?v :variant/type ?type]
       [?v :variant/gene ?g]
       [?g :gene/hugo ?gene-hugo]
       [?m :measurement/sample ?s]
       [?s :sample/id ?sample-id]]
     db "Roh2017")




;;; What are the outcomes and treatment information of patients
; who have variants in PBRM1?

; NOTE the way this is written, we only get pfs and treatment if there happens to be a timepoint that is
; tied to both a co with pfs AND a treatment-regimen
; if these pieces of info are stored on different timepoints then we won't get this info - need to rewrite

(d/q '[:find ?subject ?pfs ?treatment
       :in $ ?gene
       :where
       [?g :gene/hugo ?gene]
       [?v :variant/gene ?g]
       [?m :measurement/variant ?v]
       [?m :measurement/sample ?s]
       [?s :sample/subject ?subject]
       [?subject :subject/clinical-observations ?co]
       [?co :clinical-observation/pfs ?pfs]
       [?co :clinical-observation/timepoint ?t]
       [?t :timepoint/treatment-regimen ?tr]
       [?tr :treatment-regimen/name ?treatment]]
     db "PBRM1")





;;; What are the outcomes of patients with renal cell
; carcinoma who have variants in PBRM1, when treated with anti-PD1?

; TODO [schema] Using Nivo for now - could add an or for other compounds, or better yet, add a target to drug entity
(d/q '[:find ?subject ?pfs ?responder
       :in $ ?gene ?disease ?drug
       :where
       [?g :gene/hugo ?gene]
       [?v :variant/gene ?g]
       [?m :measurement/target ?v]
       [?m :measurement/sample ?sm]
       [?sm :sample/subject ?sb]
       [?sb :subject/uid ?subject]
       [?sb :subject/disease ?disease]
       [?co :clinical-observation/subject ?sb]
       [?co :clinical-observation/pfs ?pfs]
       [?co :clinical-observation/responder ?responder]
       [?co :clinical-observation/timepoint ?t]
       [?t :timepoint/treatment-regimen ?tr]
       [?tr :treatment-regimen/name ?drug]]
     db "PBRM1" "Renal Cell Carcinoma" "Nivolumab")



;;; Give me all the gene expression and response data for patients with NSCLC treated with anti-PD1 plus anti-CTLA4

; TODO [query] make drugs input a set/list
; Note this returns the entity ID for subjects - could also return the subject/id or subject/uid
(d/q '[:find ?subject ?pfs ?responder ?gene ?fpkm
       :in $ ?disease ?drug1 ?drug2
       :where
       [?sb :subject/disease ?disease]
       [?sb :subject/uid ?subject]
       [?sm :sample/subject ?sb]
       [?m :measurement/sample ?sm]
       [?m :measurement/fpkm ?fpkm]
       [?m :measurement/target ?g]
       [?g :gene/hugo ?gene]]
     db "Non Small Cell Lung Carcinoma" "Nivolumab" "Ipilumimab")



;;; How many treatment-naive patients have flow on matched/blood tissue in the 90 days leading up to therapy start?

; TODO [query] We should make functions that also fetch samples that roughly fit the criteria but don't have a precisely known day value
; Note that this way of querying restricts to samples that have an 'offset' value, and we probably want to not do that
; Question - does this query explode to have many many entries per subject, because of the inclusion of measurement?
; TODO [schema] Re: above question - maybe we should add a pointer directly from experimental condition to sample
(d/q '[:find [?subject ...]
       :in $ ?min-day ?max-day
       :where
       (or [?t-tumor :timepoint/type :timepoint/pre-treatment]
           [?t-tumor :timepoint/type :timepoint/pre-treatment1])
       [?t-tumor :timepoint/offset ?d-tumor]
       [(< ?d-tumor ?max-day)]
       [(> ?d-tumor ?min-day)]
       [?s-tumor :sample/timepoint ?t-tumor]
       (or [?s-tumor :sample/specimen :sample.specimen/ffpe]
           [?s-tumor :sample/specimen :sample.specimen/fresh-frozen])
       [?s-tumor :sample/subject ?subject]
       [?s-blood :sample/subject ?subject]
       (or [?s-blood :sample/specimen :sample.specimen/pbmc]
           [?s-blood :sample/specimen :sample.specimen/whole-blood])
       [?s-blood :sample/timepoint ?t-blood]
       (or [?t-blood :timepoint/type :timepoint/pre-treatment]
           [?t-blood :timepoint/type :timepoint/pre-treatment1])
       [?t-blood :timepoint/offset ?d-blood]
       [(< ?d-blood ?max-day)]
       [(> ?d-blood ?min-day)]
       [?m-blood :measurement/sample ?s-blood]
       [?e-blood :measurement-set/measurements ?m-blood]
       [?a-blood :assay/measurement-sets ?e-blood]
       [?a-blood :assay/technology :technology/flow-cytometry]
       [?m-tumor :measurement/sample ?s-tumor]
       [?e-tumor :measurement-set/measurements ?m-tumor]
       [?a-tumor :assay/measurement-sets ?e-tumor]
       [?a-tumor :assay/technology :technology/flow-cytometry]]
     db -90 0)



;;; What was the difference in the number or AEs events between pre-treated and treatment-naive patients on
; CTLA-4 + PD-1 within the first 90 days of therapy?

; assuming the query-er actually wants more info than the number of events
; TODO [query or schema] the type structure on timepoint makes this a bit difficult ... this query def doesn't capture all
; TODO [schema] still need a target entity on drug to query like this, but more completely

; fetch all subjects, AEs, and therapy line - do the count yourself
(d'q '[:find ?subject ?therapy-line ?ae
       :in $ ?days-post ?drug1 ?drug2
       :where
       [?d1 :drug/name ?drug1]
       [?dr1 :drug-regimen/drug d1]
       [?d2 :drug/name ?drug2]
       [?dr2 :drug-regimen/drug d2]
       [?tr :treatment-regimen/drug-regimens ?dr1]
       [?tr :treatment-regimen/drug-regimens ?dr2]
       [?th :therapy/treatment-regimen ?tr]
       [?subject :subject/therapies ?th]
       [?th :therapy/line ?therapy-line]
       [?t :timepoint/treatment-regimen ?tr]
       [?t :timepoint/type :timepoint/post-treatment]
       [?t :timepoint/offset ?o]
       [(<= ?o ?days-post)]
       [?co :clinical-observation/timepoint ?t]
       [?co :clinical-observation/adverse-events ?ae]
       [?co :clinical-observation/subject ?subject]]
     90 "Nivolumab" "Ipilumimab")








;;; How does response vary by gender and BMI across cancer types and checkpoint blockade therapies?

(d/q '[:find ?subject ?gender ?bmi ?responder ?disease ?drug
       :in $
       :where
       [?co :clinical-observation/subject ?subject]
       [?co :clinical-observation/bmi ?bmi]
       [?co :clinical-observation/pfs ?responder]
       [?co :clinical-observation/timepoint ?t]
       [?t :timepoint/treatment-regimen ?tr]
       [?tr :treatment-regimen/drug-regimens ?dr]
       [?dr :drug-regimen/drug ?drug]
       [?subject :subject/gender ?gender]
       [?subject :subject/disease ?disease]])




;;; Give me the IHC levels of PDL1 and the gene expression of PDL1 (CD274) from metastasis of
; head and neck cancer that have copy number variations in PDL1 (CD274)

(d/q '[:find ?ihc-percent ?gx-fpkm
       :in $ ?ihc-protein ?gx-gene ?disease ?met-bool ?cnv-gene
       :where
       [?g-cnv :gene/hugo ?cnv-gene]
       [?cnv :cnv/genes ?g-cnv]
       [?m-cnv :measurement/target ?cnv]
       [?m-cnv :measurement/sample ?sm]
       [?sm :sample/metastatic ?met-bool]
       [?sm :sample/subject ?sb]
       [?sb :subject/disease ?disease]
       [?m-ihc :measurement/sample ?sm]
       [?m-ihc :measurement/target ?p]
       [?p :protein/name ?protein]
       [?m-ihc :measurement/percentage ?ihc-percent]
       [?m-gx :measurement/sample ?sm]
       [?m-gx :measurement/target ?g-meas]
       [?g-meas :gene/hugo ?gx-gene]
       [?m-gx :measurement/fpkm ?gx-fpkm]]
     "PDL1" "CD274" "Head and Neck Carcinoma" true "CD274")




;;; How many patients have matched WES and expression data on blood/tissue before and after treatment (within x days)?

; I think the right way to answer questions like this is going to be pulling all the sample and timepoint info
; and letting the query-er follow up. Also I'm ignoring the "blood/tissue" part because I don't know what it's saying.


; TODO [query] add a query rule to deal with grouping types of timepoints
; TODO [query] this query actually requires that gx and wes be run at the same timepoint on the same sample - need to
; edit for "nearby" timepoints (within x days, as query suggests)
(d/q '[:find ?subject ?sample-before-gx ?sample-after-gx ?sample-before-wes ?sample-after-wes
       :in $ ?min-day ?max-day
       :where
       (or [?t-before :timepoint/type :timepoint/pre-treatment]
           [?t-before :timepoint/type :timepoint/pre-treatment1])
       [?sample-before-gx :sample/timepoint ?t-before]
       [?sample-before-gx :sample/subject ?subject]
       [?sample-before-wes :sample/subject ?subject]
       [?sample-before-wes :sample/timepoint ?t-before]
       [?sample-after-gx :sample/subject ?subject]
       [?sample-after-gx :sample/timepoint ?t-after]
       [?sample-after-wes :sample/timepoint ?t-after]
       [?sample-after-wes :sample/subject ?subject]
       (or [?t-after :timepoint/type :timepoint/post-treatment]
           [?t-after :timepoint/type :timepoint/post-treatment1]
           [?t-after :timepoint/type :timepoint/post-treatment2]
           [?t-after :timepoint/type :timepoint/post-treatment3]
           [?t-after :timepoint/type :timepoint/post-treatment4]
           [?t-after :timepoint/type :timepoint/post-treatment5]
           [?t-after :timepoint/type :timepoint/post-treatment6]
           [?t-after :timepoint/type :timepoint/post-treatment7]
           [?t-after :timepoint/type :timepoint/post-treatment8]
           [?t-after :timepoint/type :timepoint/post-treatment9]
           [?t-after :timepoint/type :timepoint/post-treatment10])
       [?m-after-gx :measurement/sample ?sample-after-gx]
       [?e-after-gx :measurement-set/measurements ?m-after-gx]        ; is there syntax for picking just one here?
       [?a-after-gx :assay/measurement-sets ?e-after-gx]
       (or [?a-after :assay/technology :technology/nanostring]
           [?a-after :assay/technology :technology/expression-array]
           [?a-after :assay/technology :technology/RNA-seq])
       [?m-after-wes :measurement/sample ?sample-after-wes]
       [?e-after-wes :measurement-set/measurements ?m-after-wes]
       [?a-after-wes :assay/measurement-sets ?e-after-wes]
       [?a-after-wes :assay/technology :technology/WES]
       [?m-before-gx :measurement/sample ?sample-after-gx]
       [?e-before-gx :measurement-set/measurements ?m-after-gx]        ; is there syntax for picking just one here?
       [?a-before-gx :assay/measurement-sets ?e-after-gx]
       (or [?a-before :assay/technology :technology/nanostring]
           [?a-before :assay/technology :technology/expression-array]
           [?a-before :assay/technology :technology/RNA-seq])
       [?m-before-wes :measurement/sample ?sample-before-wes]
       [?e-before-wes :measurement-set/measurements ?m-before-wes]
       [?a-before-wes :assay/measurement-sets ?e-before-wes]
       [?a-before-wes :assay/technology :technology/WES]]
     db -90 0)


;;; How many patients have paired localized and metastatic biopsy samples with sequencing?
; What is the duration between these two time points?

; TODO [schema] again probably want direct link measurement-set -> sample here
(d/q '[:find ?subject ?sample-prim ?sample-met ?t-prim ?t-met ?offset-prim ?offset-met
       :in $
       :where
       [?a-prim :assay/technology :technology/WES]
       [?a-prim :assay/measurement-sets ?e-prim]
       [?e-prim :measurement-set/measurements ?m-prim]        ; is there syntax for picking just one here?
       [?m-prim :measurement/sample ?sample-prim]
       [?sample-prim :sample/timepoint ?tp-prim]
       [?tp-prim :timepoint/type ?t-prim]
       [?tp-prim :timepoint/offset ?offset-prim]
       [?sample-prim :sample/subject ?s]
       [?sample-met :sample/subject ?s]
       [?m-met :measurement/sample ?sample-met]             ; is there a syntax for picking one here?
       [?e-prim :measurement-set/measurements ?m-met]
       [?a-prim :assay/measurement-sets ?e-prim]
       [?a-prim :assay/technology :technology/WES]
       [?sample-prim :sample/timepoint ?tp-prim]
       [?tp-prim :timepoint/type ?t-prim]
       [?tp-prim :timepoint/offset ?offset-prim]])




;;;;;;;;;;; Below are queries that we need to add things to the schema to answer


;;; How many patients who experienced toxicity on any given therapy had a microbiome sample profiled within
; 10 days of the toxicity event?
; TODO [schema] add microbiome

;;; How many with acquired resistance to PD1 have biallelic disruption to B2M, JAK1, JAK2 or IFNGR1?
; Of those, how many are focal 2-copy losses vs. LOH + Frameshift indel
; TODO [schema] add resistance (maybe. where? in clinical-observation? therapy?) and LOH

;;; For patients with a PR/CR to PD-1, how many distinct TCR clones were significantly increased
; between pre- and post- treatment?
; TODO [schema] add TCR

;;; What is the overlap in variants called with Mutect2 and Strelka2 variants across all patients with WES?
; TODO [schema] add pipeline info to measurement-set

;;; What fraction of TCR clones seen at baseline that significantly increase post-treatment are also seen 6 months
; later (for patients with at least 3 TCR times points)?
; Is this different in TCR derived from peripheral blood vs. those from tumor samples?
; TODO [schema] add TCR

;;; In MIBI/Codex images from patients who are responders to PD-1, how many tumor cells have at least one
; CD3+/CD8+ T-cell as a neighbor?
; TODO [schema] add image features

;;; What are the proportions of CD8 T-cells in post-treatment tumor tissue from patients with NSCLC treated with
; any PD-1 + combination as derived from EPIC, CIBERSORT and xCell?
; TODO [schema] add cell types to ref data, and deconvolution pipeline info to measurement-set

;;; Give me the peak-heights of all chromosomal regions from ATAC-seq that are known enhancers or promoters
; regulating CD8 T-cell function in patients who experience thyroiditis that begins 6 weeks or less after
; treatment with neoadjuvant ipilimumab for stage IIIa locally advanced metastatic melanoma.
; TODO [schema] add ATACseq

;;; Give me the intensities for all cells from CyTOF samples from pre-treatment peripheral blood, and include
; treatment information and response as well.
; TODO [schema] we don't want to add single-cell data probably. how to deal with this? ignore? pointer to raw?

;;; Were there differences in the ratio of CD8+ or CD4+ Teff : Tregs at baseline and on treatment in patients
; that had a grade3-4 irAE versus those that didn’t.
; TODO [schema] add cell types and irAE grade info (adverse-event entity in clinical?) - applies to all below

;;; Were there differences in the frequency of any immune cell subset between responders with a grade 3-4
; irAE and responders with no irAE following CTLA4+PD1 treatment

;;; Do patients that have a grade 3-4 irAE have differences in phenotype of Treg or NK cells compared to no irAE

;;; Do patients with a grade 3-4 irAE have differences in suppressor myeloid cell signatures - flow / cytof
; pbld or transcriptome tumor compared with no irAE

;;; How do myeloid, NK, cell signatures compare in treatment naive versus heavily treated patients within a
; tumor type, age and sex matched.

;;; Can you identify patients that are long term responders to antipd1 versus lose response/progress with pd1.
; What are the clinical characteristics compared to actla4 monotherapy (?don’t know if we’ll have this group;
; historical datasets) versus pd1 + ctla4 therapy.
; Comparisons in these group could include : pdl1 tumor status,
; frequency of tumor associated CD8+ PD1+ t cells at baseline and with response. Frequency of ‘exhausted’


