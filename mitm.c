#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>

#include<limits.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>

typedef uint64_t u64;       /* portable 64-bit integer */
typedef uint32_t u32;       /* portable 32-bit integer */
struct __attribute__ ((packed)) entry { u32 k; u64 v; };  /* hash table entry */

/***************************** global variables ******************************/

typedef struct {
    u64 data_val;  // La valeur hashée (f(x) ou g(y)) -> servira de CLÉ
    u64 candidate; // La pré-image (x ou y) -> servira de VALEUR
} entry_data;

typedef struct {
    entry_data *buffer; // Le tableau dynamique
    int count;          // Nombre d'éléments stockés par CE thread
    int capacity;       // Taille allouée pour CE thread
} thread_buffer_t;

u64 n = 0;         /* block size (in bits) */
u64 mask;          /* this is 2**n - 1 */

u64 dict_size;     /* number of slots in the hash table */
struct entry *A;   /* the hash table */

/* (P, C) : two plaintext-ciphertext pairs */
u32 P[2][2] = {{0, 0}, {0xffffffff, 0xffffffff}};
u32 C[2][2];

/************************ tools and utility functions *************************/

double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1E6;
}

// murmur64 hash functions, tailorized for 64-bit ints / Cf. Daniel Lemire
u64 murmur64(u64 x)
{
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdull;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ull;
    x ^= x >> 33;
    return x;
}

/* represent n in 4 bytes */
void human_format(u64 n, char *target)
{
    if (n < 1000) {
        sprintf(target, "%" PRId64, n);
        return;
    }
    if (n < 1000000) {
        sprintf(target, "%.1fK", n / 1e3);
        return;
    }
    if (n < 1000000000) {
        sprintf(target, "%.1fM", n / 1e6);
        return;
    }
    if (n < 1000000000000ll) {
        sprintf(target, "%.1fG", n / 1e9);
        return;
    }
    if (n < 1000000000000000ll) {
        sprintf(target, "%.1fT", n / 1e12);
        return;
    }
}

/******************************** SPECK block cipher **************************/

#define ROTL32(x,r) (((x)<<(r)) | (x>>(32-(r))))
#define ROTR32(x,r) (((x)>>(r)) | ((x)<<(32-(r))))

#define ER32(x,y,k) (x=ROTR32(x,8), x+=y, x^=k, y=ROTL32(y,3), y^=x)
#define DR32(x,y,k) (y^=x, y=ROTR32(y,3), x^=k, x-=y, x=ROTL32(x,8))

void Speck64128KeySchedule(const u32 K[],u32 rk[])
{
    u32 i,D=K[3],C=K[2],B=K[1],A=K[0];
    for(i=0;i<27;){
        rk[i]=A; ER32(B,A,i++);
        rk[i]=A; ER32(C,A,i++);
        rk[i]=A; ER32(D,A,i++);
    }
}

void Speck64128Encrypt(const u32 Pt[], u32 Ct[], const u32 rk[])
{
    u32 i;
    Ct[0]=Pt[0]; Ct[1]=Pt[1];
    for(i=0;i<27;)
        ER32(Ct[1],Ct[0],rk[i++]);
}

void Speck64128Decrypt(u32 Pt[], const u32 Ct[], u32 const rk[])
{
    int i;
    Pt[0]=Ct[0]; Pt[1]=Ct[1];
    for(i=26;i>=0;)
        DR32(Pt[1],Pt[0],rk[i--]);
}

/******************************** dictionary ********************************/

/*
 * "classic" hash table for 64-bit key-value pairs, with linear probing.  
 * It operates under the assumption that the keys are somewhat random 64-bit integers.
 * The keys are only stored modulo 2**32 - 5 (a prime number), and this can lead 
 * to some false positives.
 */
static const u32 EMPTY = 0xffffffff;
static const u64 PRIME = 0xfffffffb;

/* allocate a hash table with `size` slots (12*size bytes) */
void dict_setup(u64 size)
{
	dict_size = size;
	char hdsize[8];
	human_format(dict_size * sizeof(*A), hdsize);
	printf("Dictionary size: %sB\n", hdsize);

	A = malloc(sizeof(*A) * dict_size);
	if (A == NULL)
		err(1, "impossible to allocate the dictionnary");
	for (u64 i = 0; i < dict_size; i++)
		A[i].k = EMPTY;
}

void dict_setup_mpi(u64 global_size, int rank, int nprocs)
{
    // La taille locale à chaque noeud
    u64 local_size = global_size / nprocs;
    

    dict_size = local_size; 

    char hdsize[8];
    human_format(dict_size * sizeof(*A), hdsize);
    
    // Seul le rang 0 affiche l'info
    if (rank == 0) {
        printf("Total Dictionary size distributed over %d procs.\n", nprocs);
        printf("Local Dictionary size: %sB per node\n", hdsize);
    }

    A = malloc(sizeof(*A) * dict_size);
    if (A == NULL)
        err(1, "impossible to allocate the dictionnary");
    
    #pragma omp parallel for
    for (u64 i = 0; i < dict_size; i++)
        A[i].k = EMPTY;
}

/* Insert the binding key |----> value in the dictionnary */
void dict_insert(u64 key, u64 value)
{
    u64 h = murmur64(key) % dict_size;
    for (;;) {
        if (A[h].k == EMPTY)
            break;
        h += 1;
        if (h == dict_size)
            h = 0;
    }
    assert(A[h].k == EMPTY);
    A[h].k = key % PRIME;
    A[h].v = value;
}

/* Query the dictionnary with this `key`.  Write values (potentially) 
 *  matching the key in `values` and return their number. The `values`
 *  array must be preallocated of size (at least) `maxval`.
 *  The function returns -1 if there are more than `maxval` results.
 */
int dict_probe(u64 key, int maxval, u64 values[])
{
    u32 k = key % PRIME;
    u64 h = murmur64(key) % dict_size;
    int nval = 0;
    for (;;) {
        if (A[h].k == EMPTY)
            return nval;
        if (A[h].k == k) {
        	if (nval == maxval)
        		return -1;
            values[nval] = A[h].v;
            nval += 1;
        }
        h += 1;
        if (h == dict_size)
            h = 0;
   	}
}

/***************************** MITM problem ***********************************/

/* f : {0, 1}^n --> {0, 1}^n.  Speck64-128 encryption of P[0], using k */
u64 f(u64 k)
{
    assert((k & mask) == k);
    u32 K[4] = {k & 0xffffffff, k >> 32, 0, 0};
    u32 rk[27];
    Speck64128KeySchedule(K, rk);
    u32 Ct[2];
    Speck64128Encrypt(P[0], Ct, rk);
    return ((u64) Ct[0] ^ ((u64) Ct[1] << 32)) & mask;
}

/* g : {0, 1}^n --> {0, 1}^n.  speck64-128 decryption of C[0], using k */
u64 g(u64 k)
{
    assert((k & mask) == k);
    u32 K[4] = {k & 0xffffffff, k >> 32, 0, 0};
    u32 rk[27];
    Speck64128KeySchedule(K, rk);
    u32 Pt[2];
    Speck64128Decrypt(Pt, C[0], rk);
    return ((u64) Pt[0] ^ ((u64) Pt[1] << 32)) & mask;
}

bool is_good_pair(u64 k1, u64 k2)
{
    u32 Ka[4] = {k1 & 0xffffffff, k1 >> 32, 0, 0};
    u32 Kb[4] = {k2 & 0xffffffff, k2 >> 32, 0, 0};
    u32 rka[27];
    u32 rkb[27];
    Speck64128KeySchedule(Ka, rka);
    Speck64128KeySchedule(Kb, rkb);
    u32 mid[2];
    u32 Ct[2];
    Speck64128Encrypt(P[1], mid, rka);
    Speck64128Encrypt(mid, Ct, rkb);
    return (Ct[0] == C[1][0]) && (Ct[1] == C[1][1]);
}

/******************************************************************************/

/* search the "golden collision" */

/*int golden_claw_search2(int maxres, u64 k1[], u64 k2[], MPI_Comm comm)
{
    double start_time = wtime();
    u64 N = 1ull << n;

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    // --- A. Configuration OpenMP ---
    // On demande à OpenMP : "Combien de threads je peux utiliser au max ?"
    int max_threads = omp_get_max_threads();
    if (rank == 0) printf("Hybrid Mode: %d MPI procs, using %d OMP threads per proc.\n", nprocs, max_threads);

    // 1. Découpage de l'espace de recherche (Niveau MPI)
    // Chaque processus MPI s'occupe d'une plage [start_idx, end_idx]
    u64 range = N / nprocs;
    u64 remainder = N % nprocs;
    u64 start_idx = rank * range + (rank < remainder ? rank : remainder);
    u64 my_count = range + (rank < remainder ? 1 : 0);
    u64 end_idx = start_idx + my_count;

    // --- PHASE 1 : REMPLISSAGE (FILL) ---
    
    // Allocation du tableau principal qui contiendra les buffers de TOUS les threads.
    // thread_buffers[id_thread] donnera accès aux buffers du thread numéro 'id_thread'.
    thread_buffer_t **thread_buffers = malloc(max_threads * sizeof(thread_buffer_t*));
    
    // --- DÉBUT DE LA ZONE PARALLÈLE ---
    // À partir d'ici, le code est exécuté par 'max_threads' threads en même temps.
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();     // "Qui suis-je ?" (0, 1, 2...)
        int n_threads = omp_get_num_threads(); // "Combien sommes-nous réellement ?"

        // 2. Allocation locale : Chaque thread prépare ses propres "bacs"
        // On alloue une ligne pour ce thread, contenant 'nprocs' buffers (un par destination MPI)
        thread_buffers[tid] = malloc(nprocs * sizeof(thread_buffer_t));
        
        // Petite estimation mathématique pour la taille initiale des buffers
        // (Total à traiter) / (Nb destinations) / (Nb threads) + marge de sécurité
        int est_size = (my_count / nprocs / n_threads) * 1.5 + 128;

        for (int i = 0; i < nprocs; i++) {
            thread_buffers[tid][i].capacity = est_size;
            thread_buffers[tid][i].count = 0;
            thread_buffers[tid][i].buffer = malloc(est_size * sizeof(entry_data));
            // (Idéalement, ajouter une vérification si malloc échoue ici)
        }

        // 3. La Boucle de Remplissage Parallèle
        // #pragma omp for : OpenMP divise automatiquement la boucle entre les threads.
        // schedule(static) : Donne des blocs contigus de 'x' à chaque thread (meilleur pour le cache CPU)
        #pragma omp for schedule(static)
        for (u64 x = start_idx; x < end_idx; x++) {
            u64 z = f(x);               // Calcul lourd
            u64 hash = murmur64(z);
            int target = hash % nprocs; // Destination MPI
            
            // Accès au buffer PRIVÉ du thread 'tid' pour la destination 'target'
            thread_buffer_t *tb = &thread_buffers[tid][target];

            // Gestion de la capacité (realloc si plein)
            // Pas besoin de verrous ici car seul CE thread touche à CE buffer !
            if (tb->count >= tb->capacity) {
                tb->capacity *= 2;
                tb->buffer = realloc(tb->buffer, tb->capacity * sizeof(entry_data));
            }

            // Stockage
            tb->buffer[tb->count].data_val = z;
            tb->buffer[tb->count].candidate = x;
            tb->count++;
        }
    }
    // --- FIN DE LA ZONE PARALLÈLE ---
    // Ici, une "barrière implicite" attend que tous les threads aient fini.
    
    // ... LA SUITE (Fusion des buffers) VIENDRA À L'ÉTAPE 3 ...

    // --- SUITE DE L'ÉTAPE 2 (Juste après l'accolade fermante du pragma parallel) ---

    // --- ÉTAPE 3 : FUSION DES DONNÉES (MERGE) ---
    // On doit rassembler les données de tous les threads pour MPI.

    // 1. Calculer combien d'éléments au total on doit envoyer à chaque processus MPI
    // On additionne les compteurs de tous les threads pour chaque destination.
    int *send_counts = calloc(nprocs, sizeof(int));
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            send_counts[p] += thread_buffers[t][p].count;
        }
    }

    // 2. Calculer les positions (offsets) dans le futur grand buffer unique
    // send_displs[i] indique où commencent les données pour le processus 'i'
    int *send_displs = malloc(nprocs * sizeof(int));
    u64 total_send = 0;
    
    send_displs[0] = 0;
    total_send += send_counts[0];
    
    for (int i = 1; i < nprocs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
        total_send += send_counts[i];
    }

    // 3. Allouer le grand buffer unique ("Flat Buffer") qui sera envoyé via MPI
    entry_data *flat_send_buf = malloc(total_send * sizeof(entry_data));
    
    // 4. Copie effective : On vide les sacs des threads dans le grand container
    // 'current_offsets' nous aide à savoir où écrire pour chaque destination
    int *current_offsets = calloc(nprocs, sizeof(int));
    
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            int count = thread_buffers[t][p].count;
            if (count > 0) {
                // Calcul de la position d'écriture : Début de la zone P + Décalage actuel
                u64 dest_idx = send_displs[p] + current_offsets[p];
                
                // Copie rapide (memcpy est optimisé)
                memcpy(&flat_send_buf[dest_idx], thread_buffers[t][p].buffer, count * sizeof(entry_data));
                
                // On met à jour le décalage pour le prochain thread qui voudra écrire pour P
                current_offsets[p] += count;
                
                // IMPORTANT : On libère la mémoire du thread dès qu'on a copié !
                // Cela évite d'avoir les données en double en RAM trop longtemps.
                free(thread_buffers[t][p].buffer); 
            }
        }
        free(thread_buffers[t]); // On libère la structure du thread t
    }
    free(thread_buffers);  // On libère le tableau de pointeurs
    free(current_offsets); // On libère l'outil temporaire

    // --- FIN DE L'ÉTAPE 3 ---

    double mid_time = wtime();
    // Petit print de debug (seulement le rank 0 pour pas spammer)
    if (rank == 0) printf("Fill preparation done in %.2fs. Ready to exchange.\n", mid_time - start_time);

    
    // --- ÉTAPE 3.5 : ÉCHANGE MPI ET INSERTION (Code MPI standard) ---
    
    // 1. Échanger les tailles (Combien vais-je recevoir ?)
    int *recv_counts = malloc(nprocs * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    // 2. Calculer les déplacements de réception
    int *recv_displs = malloc(nprocs * sizeof(int));
    u64 total_recv = 0;
    recv_displs[0] = 0;
    total_recv = recv_counts[0];
    
    for (int i = 1; i < nprocs; i++) {
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        total_recv += recv_counts[i];
    }

    // 3.4 Allouer le buffer de réception
    entry_data *recv_buf = malloc(total_recv * sizeof(entry_data));
    if (recv_buf == NULL && total_recv > 0) {
         err(1, "Malloc recv_buf failed");
    }

    // 3.5 Conversion en BYTES pour l'envoi
    // MPI gère mal les types custom si on ne crée pas un MPI_Datatype.
    // Astuce classique : on triche en disant qu'on envoie des MPI_BYTE (octets).
    // Il faut donc multiplier tous les counts et displs par sizeof(entry_data).
    
    int *s_cnts_b = malloc(nprocs * sizeof(int));
    int *r_cnts_b = malloc(nprocs * sizeof(int));
    int *s_displs_b = malloc(nprocs * sizeof(int));
    int *r_displs_b = malloc(nprocs * sizeof(int));

    for(int i=0; i<nprocs; i++) {
        s_cnts_b[i] = send_counts[i] * sizeof(entry_data);
        r_cnts_b[i] = recv_counts[i] * sizeof(entry_data);
        s_displs_b[i] = send_displs[i] * sizeof(entry_data);
        r_displs_b[i] = recv_displs[i] * sizeof(entry_data);
    }

    // 3.6 L'échange MPI (enfin !)
    MPI_Alltoallv(flat_send_buf, s_cnts_b, s_displs_b, MPI_BYTE,
                  recv_buf,      r_cnts_b, r_displs_b, MPI_BYTE,
                  comm);

    // 3.7 Insertion dans le Dictionnaire Local
    // Maintenant, 'recv_buf' contient toutes les données dont JE suis responsable.
    
    for (u64 i = 0; i < total_recv; i++) {
        // recv_buf[i].data_val  c'est le hash z (la CLÉ)
        // recv_buf[i].candidate c'est la valeur x
        dict_insert(recv_buf[i].data_val, recv_buf[i].candidate);
    }
    
    // 3.8 Nettoyage mémoire de cette phase
    free(flat_send_buf);
    free(recv_buf);
    free(send_counts);
    free(recv_counts);
    free(send_displs);
    free(recv_displs);

    // On garde les buffers _b pour la phase suivante ou on les free et realloc ? 
    // Mieux vaut tout free pour être propre, on réallouera pour la phase Probe.
    free(s_cnts_b); free(r_cnts_b); free(s_displs_b); free(r_displs_b);

    MPI_Barrier(comm); // Attendre que tout le monde ait fini de remplir
    mid_time = wtime();
    if (rank == 0) printf("Fill phase logic done in %.1fs\n", mid_time - start_time);


    // --- PHASE 2 : SONDAGE (PROBE) ---
    // Réallocation des structures pour la phase 2 (on repart à zéro)
    
    // 1. Réallocation de la structure pour les threads (car libérée en Phase 1)
    thread_buffers = malloc(max_threads * sizeof(thread_buffer_t*));
    
    // 2. Zone Parallèle : Remplissage des buffers avec g(y)
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int n_threads = omp_get_num_threads();

        // Chaque thread alloue sa ligne de buffers
        thread_buffers[tid] = malloc(nprocs * sizeof(thread_buffer_t));
        
        // Estimation de la taille nécessaire
        int est_size = (my_count / nprocs / n_threads) * 1.5 + 128;

        for (int i = 0; i < nprocs; i++) {
            thread_buffers[tid][i].capacity = est_size;
            thread_buffers[tid][i].count = 0;
            thread_buffers[tid][i].buffer = malloc(est_size * sizeof(entry_data));
        }

        // Boucle Parallèle sur les Y
        // On utilise schedule(static) pour découper la plage [start_idx, end_idx] équitablement
        #pragma omp for schedule(static)
        for (u64 y_val = start_idx; y_val < end_idx; y_val++) {
            u64 z_prime = g(y_val);     // Calcul de g(y)
            u64 hash = murmur64(z_prime);
            int target = hash % nprocs; // Même logique de hachage que Phase 1
            
            thread_buffer_t *tb = &thread_buffers[tid][target];

            // Gestion dynamique de la mémoire (sans conflit car buffer privé)
            if (tb->count >= tb->capacity) {
                tb->capacity *= 2;
                tb->buffer = realloc(tb->buffer, tb->capacity * sizeof(entry_data));
            }

            // On stocke la paire (g(y), y)
            tb->buffer[tb->count].data_val = z_prime;
            tb->buffer[tb->count].candidate = y_val;
            tb->count++;
        }
    } // Fin de la zone parallèle

    // 3. Fusion (Merge) : On rassemble tout pour MPI
    
    // Remise à zéro des compteurs (Important !)
    //memset(send_counts, 0, nprocs * sizeof(int));
    send_counts = calloc(nprocs, sizeof(int));

    // Calcul des totaux par destination
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            send_counts[p] += thread_buffers[t][p].count;
        }
    }

    // Recalcul des offsets pour le grand buffer
    send_displs = malloc(nprocs * sizeof(int));
    total_send = 0;
    send_displs[0] = 0;
    total_send += send_counts[0];
    for (int i = 1; i < nprocs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
        total_send += send_counts[i];
    }

    // Réallocation du buffer plat
    flat_send_buf = malloc(total_send * sizeof(entry_data));
    
    // Reset des offsets temporaires
    //memset(current_offsets, 0, nprocs * sizeof(int));
    current_offsets = calloc(nprocs, sizeof(int));

    // Copie des données
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            int count = thread_buffers[t][p].count;
            if (count > 0) {
                u64 dest_idx = send_displs[p] + current_offsets[p];
                memcpy(&flat_send_buf[dest_idx], thread_buffers[t][p].buffer, count * sizeof(entry_data));
                current_offsets[p] += count;
                free(thread_buffers[t][p].buffer); // Libération immédiate
            }
        }
        free(thread_buffers[t]);
    }
    free(thread_buffers);
    free(current_offsets);

    // 4. Échange MPI (Probe)

    recv_counts = malloc(nprocs * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    // Recalcul recv_displs
    recv_displs = malloc(nprocs * sizeof(int));
    total_recv = 0;
    recv_displs[0] = 0;
    total_recv = recv_counts[0];
    for (int i = 1; i < nprocs; i++) {
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        total_recv += recv_counts[i];
    }

    // Allocation recv_buf
    recv_buf = malloc(total_recv * sizeof(entry_data));

    // Mise à jour des structures MPI_BYTE
    s_cnts_b = malloc(nprocs * sizeof(int));
    r_cnts_b = malloc(nprocs * sizeof(int));
    s_displs_b = malloc(nprocs * sizeof(int));
    r_displs_b = malloc(nprocs * sizeof(int));

    for(int i=0; i<nprocs; i++) {
        s_cnts_b[i] = send_counts[i] * sizeof(entry_data);
        r_cnts_b[i] = recv_counts[i] * sizeof(entry_data);
        s_displs_b[i] = send_displs[i] * sizeof(entry_data);
        r_displs_b[i] = recv_displs[i] * sizeof(entry_data);
    }

    // L'échange
    MPI_Alltoallv(flat_send_buf, s_cnts_b, s_displs_b, MPI_BYTE,
                  recv_buf,      r_cnts_b, r_displs_b, MPI_BYTE,
                  comm);

    // 5. Vérification Parallèle des Résultats (Probe Final)
    
    int nres = 0;
    
    #pragma omp parallel 
    {
        // Buffer déclaré ICI pour être privé à chaque thread
        u64 x_candidates[256]; 
        
        #pragma omp for
        for (u64 i = 0; i < total_recv; i++) {
            u64 val_to_find = recv_buf[i].data_val;  // g(y)
            u64 candidate_y = recv_buf[i].candidate; // y

            // Lecture parallèle (Thread-Safe)
            int nx = dict_probe(val_to_find, 256, x_candidates);
            
            for (int k = 0; k < nx; k++) {
                if (is_good_pair(x_candidates[k], candidate_y)) {
                    
                    // SECTION CRITIQUE : Un seul thread écrit le résultat
                    #pragma omp critical
                    {
                        if (nres < maxres) {
                            k1[nres] = x_candidates[k];
                            k2[nres] = candidate_y;
                            nres++;
                            printf("#####################################################\n");
                            printf("[Rank %d Thread %d] FOUND SOLUTION: k1=%" PRIx64 ", k2=%" PRIx64 "\n", 
                                    rank, omp_get_thread_num(), x_candidates[k], candidate_y);
                            printf("#####################################################\n");
                        }
                    }
                }
            }
        }
    }

    // --- Nettoyage Final ---
    free(flat_send_buf); free(recv_buf);
    free(send_counts); free(recv_counts);
    free(send_displs); free(recv_displs);
    // free(send_capacities); // Plus utilisé ici car géré dans les threads
    free(s_cnts_b); free(r_cnts_b); free(s_displs_b); free(r_displs_b);

    return nres;
}*/

int golden_claw_search2(int maxres, u64 k1[], u64 k2[], MPI_Comm comm)
{
    double start_time = wtime();
    u64 N = 1ull << n;

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int max_threads = omp_get_max_threads();
    if (rank == 0) printf("Hybrid Mode: %d MPI procs, using %d OMP threads per proc.\n", nprocs, max_threads);

    
    MPI_Datatype dt_entry;
    MPI_Type_contiguous(sizeof(entry_data), MPI_BYTE, &dt_entry);
    MPI_Type_commit(&dt_entry);

    u64 range = N / nprocs;
    u64 remainder = N % nprocs;
    u64 start_idx = rank * range + (rank < remainder ? rank : remainder);
    u64 my_count = range + (rank < remainder ? 1 : 0);
    u64 end_idx = start_idx + my_count;

    // PHASE 1 : REMPLISSAGE (FILL)
    thread_buffer_t **thread_buffers = malloc(max_threads * sizeof(thread_buffer_t*));
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int n_threads = omp_get_num_threads();
        thread_buffers[tid] = malloc(nprocs * sizeof(thread_buffer_t));
        int est_size = (my_count / nprocs / n_threads) * 1.5 + 128;

        for (int i = 0; i < nprocs; i++) {
            thread_buffers[tid][i].capacity = est_size;
            thread_buffers[tid][i].count = 0;
            thread_buffers[tid][i].buffer = malloc(est_size * sizeof(entry_data));
        }

        #pragma omp for schedule(static)
        for (u64 x = start_idx; x < end_idx; x++) {
            u64 z = f(x);
            u64 hash = murmur64(z);
            int target = hash % nprocs;
            
            thread_buffer_t *tb = &thread_buffers[tid][target];
            if (tb->count >= tb->capacity) {
                tb->capacity *= 2;
                tb->buffer = realloc(tb->buffer, tb->capacity * sizeof(entry_data));
            }
            tb->buffer[tb->count].data_val = z;
            tb->buffer[tb->count].candidate = x;
            tb->count++;
        }
    } 

    // Fusion Phase 1
    int *send_counts = calloc(nprocs, sizeof(int));
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            send_counts[p] += thread_buffers[t][p].count;
        }
    }

    int *send_displs = malloc(nprocs * sizeof(int));
    u64 total_send = 0;
    send_displs[0] = 0;
    total_send += send_counts[0];
    for (int i = 1; i < nprocs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
        total_send += send_counts[i];
    }

    entry_data *flat_send_buf = malloc(total_send * sizeof(entry_data));
    int *current_offsets = calloc(nprocs, sizeof(int));
    
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            int count = thread_buffers[t][p].count;
            if (count > 0) {
                u64 dest_idx = send_displs[p] + current_offsets[p];
                memcpy(&flat_send_buf[dest_idx], thread_buffers[t][p].buffer, count * sizeof(entry_data));
                current_offsets[p] += count;
                free(thread_buffers[t][p].buffer); 
            }
        }
        free(thread_buffers[t]);
    }
    free(thread_buffers);
    free(current_offsets);

    double mid_time = wtime();
    // Temps Fill
    if (rank == 0) printf("Fill: %.1fs\n", mid_time - start_time); 

    // Échange MPI Phase 1
    int *recv_counts = malloc(nprocs * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    int *recv_displs = malloc(nprocs * sizeof(int));
    u64 total_recv = 0;
    recv_displs[0] = 0;
    total_recv = recv_counts[0];
    for (int i = 1; i < nprocs; i++) {
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        total_recv += recv_counts[i];
    }

    entry_data *recv_buf = malloc(total_recv * sizeof(entry_data));
    MPI_Alltoallv(flat_send_buf, send_counts, send_displs, dt_entry,
                  recv_buf,      recv_counts, recv_displs, dt_entry, comm);

    for (u64 i = 0; i < total_recv; i++) {
        dict_insert(recv_buf[i].data_val, recv_buf[i].candidate);
    }

    free(flat_send_buf); free(recv_buf);
    free(send_counts); free(recv_counts);
    free(send_displs); free(recv_displs);

    MPI_Barrier(comm);

    //PHASE 2 : SONDAGE (PROBE)
    
    thread_buffers = malloc(max_threads * sizeof(thread_buffer_t*));
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int n_threads = omp_get_num_threads();
        thread_buffers[tid] = malloc(nprocs * sizeof(thread_buffer_t));
        int est_size = (my_count / nprocs / n_threads) * 1.5 + 128;

        for (int i = 0; i < nprocs; i++) {
            thread_buffers[tid][i].capacity = est_size;
            thread_buffers[tid][i].count = 0;
            thread_buffers[tid][i].buffer = malloc(est_size * sizeof(entry_data));
        }

        #pragma omp for schedule(static)
        for (u64 y_val = start_idx; y_val < end_idx; y_val++) {
            u64 z_prime = g(y_val); 
            u64 hash = murmur64(z_prime);
            int target = hash % nprocs;
            
            thread_buffer_t *tb = &thread_buffers[tid][target];
            if (tb->count >= tb->capacity) {
                tb->capacity *= 2;
                tb->buffer = realloc(tb->buffer, tb->capacity * sizeof(entry_data));
            }
            tb->buffer[tb->count].data_val = z_prime;
            tb->buffer[tb->count].candidate = y_val;
            tb->count++;
        }
    }

    // Fusion Phase 2
    send_counts = calloc(nprocs, sizeof(int));
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            send_counts[p] += thread_buffers[t][p].count;
        }
    }

    send_displs = malloc(nprocs * sizeof(int));
    total_send = 0;
    send_displs[0] = 0;
    total_send += send_counts[0];
    for (int i = 1; i < nprocs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
        total_send += send_counts[i];
    }

    flat_send_buf = malloc(total_send * sizeof(entry_data));
    current_offsets = calloc(nprocs, sizeof(int));
    
    for (int t = 0; t < max_threads; t++) {
        for (int p = 0; p < nprocs; p++) {
            int count = thread_buffers[t][p].count;
            if (count > 0) {
                u64 dest_idx = send_displs[p] + current_offsets[p];
                memcpy(&flat_send_buf[dest_idx], thread_buffers[t][p].buffer, count * sizeof(entry_data));
                current_offsets[p] += count;
                free(thread_buffers[t][p].buffer);
            }
        }
        free(thread_buffers[t]);
    }
    free(thread_buffers);
    free(current_offsets);

    // Échange MPI Phase 2
    recv_counts = malloc(nprocs * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    recv_displs = malloc(nprocs * sizeof(int));
    total_recv = 0;
    recv_displs[0] = 0;
    total_recv = recv_counts[0];
    for (int i = 1; i < nprocs; i++) {
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        total_recv += recv_counts[i];
    }

    recv_buf = malloc(total_recv * sizeof(entry_data));
    MPI_Alltoallv(flat_send_buf, send_counts, send_displs, dt_entry,
                  recv_buf,      recv_counts, recv_displs, dt_entry, comm);

    //VÉRIFICATION PARALLÈLE
    int nres = 0;
    u64 ncandidates = 0;
    
    
    #pragma omp parallel reduction(+:ncandidates)
    {
        u64 x_candidates[256]; 
        
        #pragma omp for
        for (u64 i = 0; i < total_recv; i++) {
            u64 val_to_find = recv_buf[i].data_val;
            u64 candidate_y = recv_buf[i].candidate;

            int nx = dict_probe(val_to_find, 256, x_candidates);
            
            ncandidates += nx;
            
            for (int k = 0; k < nx; k++) {
                if (is_good_pair(x_candidates[k], candidate_y)) {
                    #pragma omp critical
                    {
                        if (nres < maxres) {
                            k1[nres] = x_candidates[k];
                            k2[nres] = candidate_y;
                            nres++;
                            printf("#####################################################\n");
                            printf("[Rank %d Thread %d] FOUND SOLUTION: k1=%" PRIx64 ", k2=%" PRIx64 "\n", 
                                    rank, omp_get_thread_num(), x_candidates[k], candidate_y);
                            printf("#####################################################\n");
                        }
                    }
                }
            }
        }
    }

    // On somme les ncandidates de tous les processus pour avoir le total global
    u64 total_candidates_global = 0;
    MPI_Reduce(&ncandidates, &total_candidates_global, 1, MPI_UINT64_T, MPI_SUM, 0, comm);

    if (rank == 0) {
        printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid_time, total_candidates_global);
    }

    // Nettoyage Final
    MPI_Type_free(&dt_entry);
    free(flat_send_buf); free(recv_buf);
    free(send_counts); free(recv_counts);
    free(send_displs); free(recv_displs);

    return nres;
}


int golden_claw_search3(int maxres, u64 k1[], u64 k2[], MPI_Comm comm)
{
    double start_time = wtime();
    u64 N = 1ull << n;

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    int max_threads = omp_get_max_threads();
    if (rank == 0) printf("Chunking Mode: %d MPI procs, %d OMP threads.\n", nprocs, max_threads);

    MPI_Datatype dt_entry;
    MPI_Type_contiguous(sizeof(entry_data), MPI_BYTE, &dt_entry);
    MPI_Type_commit(&dt_entry);

    u64 range = N / nprocs;
    u64 remainder = N % nprocs;
    u64 global_start_idx = rank * range + (rank < remainder ? rank : remainder);
    u64 my_total_count = range + (rank < remainder ? 1 : 0);
    
    int NUM_CHUNKS = 256; 
    u64 chunk_size = my_total_count / NUM_CHUNKS;

    //PHASE 1 : REMPLISSAGE (FILL)

    for (int c = 0; c < NUM_CHUNKS; c++) {
        
        u64 c_start = global_start_idx + c * chunk_size;
        u64 c_count = (c == NUM_CHUNKS - 1) ? (my_total_count - c * chunk_size) : chunk_size;
        u64 c_end = c_start + c_count;

        thread_buffer_t **thread_buffers = malloc(max_threads * sizeof(thread_buffer_t*));
        
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int n_threads = omp_get_num_threads();
            thread_buffers[tid] = malloc(nprocs * sizeof(thread_buffer_t));
            
            int est_size = (c_count / nprocs / n_threads) * 1.5 + 128;

            for (int i = 0; i < nprocs; i++) {
                thread_buffers[tid][i].capacity = est_size;
                thread_buffers[tid][i].count = 0;
                thread_buffers[tid][i].buffer = malloc(est_size * sizeof(entry_data));
            }

            #pragma omp for schedule(static)
            for (u64 x = c_start; x < c_end; x++) {
                u64 z = f(x);
                u64 hash = murmur64(z);
                int target = hash % nprocs;
                
                thread_buffer_t *tb = &thread_buffers[tid][target];
                if (tb->count >= tb->capacity) {
                    tb->capacity *= 2;
                    tb->buffer = realloc(tb->buffer, tb->capacity * sizeof(entry_data));
                }
                tb->buffer[tb->count].data_val = z;
                tb->buffer[tb->count].candidate = x;
                tb->count++;
            }
        } 

        // Fusion Chunk
        int *send_counts = calloc(nprocs, sizeof(int));
        for (int t = 0; t < max_threads; t++) {
            for (int p = 0; p < nprocs; p++) {
                send_counts[p] += thread_buffers[t][p].count;
            }
        }
        for (int p = 0; p < nprocs; p++) {
            assert(send_counts[p] < INT_MAX);
        }

        int *send_displs = malloc(nprocs * sizeof(int));
        u64 total_send = 0;
        send_displs[0] = 0;
        total_send += send_counts[0];
        for (int i = 1; i < nprocs; i++) {
            send_displs[i] = send_displs[i-1] + send_counts[i-1];
            total_send += send_counts[i];
        }

        entry_data *flat_send_buf = malloc(total_send * sizeof(entry_data));
        int *current_offsets = calloc(nprocs, sizeof(int));
        
        for (int t = 0; t < max_threads; t++) {
            for (int p = 0; p < nprocs; p++) {
                int count = thread_buffers[t][p].count;
                if (count > 0) {
                    u64 dest_idx = send_displs[p] + current_offsets[p];
                    memcpy(&flat_send_buf[dest_idx], thread_buffers[t][p].buffer, count * sizeof(entry_data));
                    current_offsets[p] += count;
                    free(thread_buffers[t][p].buffer); 
                }
            }
            free(thread_buffers[t]);
        }
        free(thread_buffers);
        free(current_offsets);

        // Échange MPI Chunk
        int *recv_counts = malloc(nprocs * sizeof(int));
        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

        int *recv_displs = malloc(nprocs * sizeof(int));
        u64 total_recv = 0;
        recv_displs[0] = 0;
        total_recv = recv_counts[0];
        for (int i = 1; i < nprocs; i++) {
            recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            total_recv += recv_counts[i];
        }

        entry_data *recv_buf = malloc(total_recv * sizeof(entry_data));
        MPI_Alltoallv(flat_send_buf, send_counts, send_displs, dt_entry,
                      recv_buf,      recv_counts, recv_displs, dt_entry, comm);

        // Insertion Chunk
        for (u64 i = 0; i < total_recv; i++) {
            dict_insert(recv_buf[i].data_val, recv_buf[i].candidate);
        }

        // Nettoyage Chunk
        free(flat_send_buf); free(recv_buf);
        free(send_counts); free(recv_counts);
        free(send_displs); free(recv_displs);
        
        // Print pour voir l'avancement (on peut supp)
        if (rank == 0) printf("  Phase 1: Chunk %d/%d processed.\n", c+1, NUM_CHUNKS);
    }

    double mid_time = wtime();
    if (rank == 0) printf("Fill: %.1fs\n", mid_time - start_time); 

    MPI_Barrier(comm);

    //PHASE 2 : SONDAGE (PROBE)
    
    int nres = 0;
    u64 ncandidates = 0;

    for (int c = 0; c < NUM_CHUNKS; c++) {
        
        u64 c_start = global_start_idx + c * chunk_size;
        u64 c_count = (c == NUM_CHUNKS - 1) ? (my_total_count - c * chunk_size) : chunk_size;
        u64 c_end = c_start + c_count;

        thread_buffer_t **thread_buffers = malloc(max_threads * sizeof(thread_buffer_t*));
        
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int n_threads = omp_get_num_threads();
            thread_buffers[tid] = malloc(nprocs * sizeof(thread_buffer_t));
            int est_size = (c_count / nprocs / n_threads) * 1.5 + 128;

            for (int i = 0; i < nprocs; i++) {
                thread_buffers[tid][i].capacity = est_size;
                thread_buffers[tid][i].count = 0;
                thread_buffers[tid][i].buffer = malloc(est_size * sizeof(entry_data));
            }

            #pragma omp for schedule(static)
            for (u64 y_val = c_start; y_val < c_end; y_val++) {
                u64 z_prime = g(y_val); 
                u64 hash = murmur64(z_prime);
                int target = hash % nprocs;
                
                thread_buffer_t *tb = &thread_buffers[tid][target];
                if (tb->count >= tb->capacity) {
                    tb->capacity *= 2;
                    tb->buffer = realloc(tb->buffer, tb->capacity * sizeof(entry_data));
                }
                tb->buffer[tb->count].data_val = z_prime;
                tb->buffer[tb->count].candidate = y_val;
                tb->count++;
            }
        }

        // Fusion Chunk Phase 2
        int *send_counts = calloc(nprocs, sizeof(int));
        for (int t = 0; t < max_threads; t++) {
            for (int p = 0; p < nprocs; p++) {
                send_counts[p] += thread_buffers[t][p].count;
            }
        }

        for (int p = 0; p < nprocs; p++) {
            assert(send_counts[p] < INT_MAX);
        }

        int *send_displs = malloc(nprocs * sizeof(int));
        u64 total_send = 0;
        send_displs[0] = 0;
        total_send += send_counts[0];
        for (int i = 1; i < nprocs; i++) {
            send_displs[i] = send_displs[i-1] + send_counts[i-1];
            total_send += send_counts[i];
        }

        entry_data *flat_send_buf = malloc(total_send * sizeof(entry_data));
        int *current_offsets = calloc(nprocs, sizeof(int));
        
        for (int t = 0; t < max_threads; t++) {
            for (int p = 0; p < nprocs; p++) {
                int count = thread_buffers[t][p].count;
                if (count > 0) {
                    u64 dest_idx = send_displs[p] + current_offsets[p];
                    memcpy(&flat_send_buf[dest_idx], thread_buffers[t][p].buffer, count * sizeof(entry_data));
                    current_offsets[p] += count;
                    free(thread_buffers[t][p].buffer);
                }
            }
            free(thread_buffers[t]);
        }
        free(thread_buffers);
        free(current_offsets);

        // Échange MPI Chunk Phase 2
        int *recv_counts = malloc(nprocs * sizeof(int));
        MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

        int *recv_displs = malloc(nprocs * sizeof(int));
        u64 total_recv = 0;
        recv_displs[0] = 0;
        total_recv = recv_counts[0];
        for (int i = 1; i < nprocs; i++) {
            recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
            total_recv += recv_counts[i];
        }

        entry_data *recv_buf = malloc(total_recv * sizeof(entry_data));
        MPI_Alltoallv(flat_send_buf, send_counts, send_displs, dt_entry,
                      recv_buf,      recv_counts, recv_displs, dt_entry, comm);

        //VÉRIFICATION CHUNK
        #pragma omp parallel reduction(+:ncandidates)
        {
            u64 x_candidates[256]; 
            
            #pragma omp for
            for (u64 i = 0; i < total_recv; i++) {
                u64 val_to_find = recv_buf[i].data_val;
                u64 candidate_y = recv_buf[i].candidate;

                int nx = dict_probe(val_to_find, 256, x_candidates);
                ncandidates += nx;
                
                for (int k = 0; k < nx; k++) {
                    if (is_good_pair(x_candidates[k], candidate_y)) {
                        #pragma omp critical
                        {
                            if (nres < maxres) {
                                k1[nres] = x_candidates[k];
                                k2[nres] = candidate_y;
                                nres++;
                                printf("#####################################################\n");
                                printf("[Rank %d Chunk %d] FOUND SOLUTION: k1=%" PRIx64 ", k2=%" PRIx64 "\n", 
                                        rank, c, x_candidates[k], candidate_y);
                                printf("#####################################################\n");
                            }
                        }
                    }
                }
            }
        }

        // Nettoyage Chunk Phase 2
        free(flat_send_buf); free(recv_buf);
        free(send_counts); free(recv_counts);
        free(send_displs); free(recv_displs);
        if (rank == 0) printf("  [Phase 2] Chunk %d/%d processed.\n", c+1, NUM_CHUNKS);

    }

    u64 total_candidates_global = 0;
    MPI_Reduce(&ncandidates, &total_candidates_global, 1, MPI_UINT64_T, MPI_SUM, 0, comm);

    if (rank == 0) {
        printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid_time, total_candidates_global);
    }

    MPI_Type_free(&dt_entry);
    return nres;
}

/************************** command-line options ****************************/

void usage(char **argv)
{
        printf("%s [OPTIONS]\n\n", argv[0]);
        printf("Options:\n");
        printf("--n N                       block size [default 24]\n");
        printf("--C0 N                      1st ciphertext (in hex)\n");
        printf("--C1 N                      2nd ciphertext (in hex)\n");
        printf("\n");
        printf("All arguments are required\n");
        exit(0);
}

void process_command_line_options(int argc, char ** argv)
{
        struct option longopts[4] = {
                {"n", required_argument, NULL, 'n'},
                {"C0", required_argument, NULL, '0'},
                {"C1", required_argument, NULL, '1'},
                {NULL, 0, NULL, 0}
        };
        char ch;
        int set = 0;
        while ((ch = getopt_long(argc, argv, "", longopts, NULL)) != -1) {
                switch (ch) {
                case 'n':
                        n = atoi(optarg);
                        mask = (1ull << n) - 1;
                        break;
                case '0':
                        set |= 1;
                        u64 c0 = strtoull(optarg, NULL, 16);
                        C[0][0] = c0 & 0xffffffff;
                        C[0][1] = c0 >> 32;
                        break;
                case '1':
                        set |= 2;
                        u64 c1 = strtoull(optarg, NULL, 16);
                        C[1][0] = c1 & 0xffffffff;
                        C[1][1] = c1 >> 32;
                        break;
                default:
                        errx(1, "Unknown option\n");
                }
        }
        if (n == 0 || set != 3) {
        	usage(argv);
        	exit(1);
        }
}

/******************************************************************************/

int main(int argc, char **argv)
{
    // Initialisation MPI
    MPI_Init(&argc, &argv);
    
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    process_command_line_options(argc, argv);

    if (rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n", 
            (int) n, C[0][0], C[0][1], C[1][0], C[1][1]);
    }

    // Allocation distribuée
    dict_setup_mpi(1.125 * (1ull << n), rank, nprocs);

    /* search */
    u64 k1[16], k2[16];
    int nkey = golden_claw_search3(16, k1, k2, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);

    /* validation */
    for (int i = 0; i < nkey; i++) {
        assert(f(k1[i]) == g(k2[i]));
        assert(is_good_pair(k1[i], k2[i]));     
        
        printf("#################################################################\n");
        printf("[Rank %d] FINAL VALIDATED SOLUTION: (%" PRIx64 ", %" PRIx64 ")\n", rank, k1[i], k2[i]);
        printf("#################################################################\n");
    }

    free(A);
    MPI_Finalize();
    
    return 0;
}