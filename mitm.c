#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <assert.h>
#include <getopt.h>
#include <err.h>
#include <assert.h>

#include <string.h>
#include <mpi.h>

typedef uint64_t u64;       /* portable 64-bit integer */
typedef uint32_t u32;       /* portable 32-bit integer */
struct __attribute__ ((packed)) entry { u32 k; u64 v; };  /* hash table entry */

/***************************** global variables ******************************/

typedef struct {
    u64 data_val;  // La valeur hashée (f(x) ou g(y)) -> servira de CLÉ
    u64 candidate; // La pré-image (x ou y) -> servira de VALEUR
} entry_data;

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

/* allocate a hash table logic for MPI */
void dict_setup_mpi(u64 global_size, int rank, int nprocs)
{
    // Chaque nœud prend sa part du gâteau
    u64 local_size = global_size / nprocs;
    
    // On ajoute un peu de marge pour gérer les arrondis et éviter les collisions
    dict_size = local_size; 

    char hdsize[8];
    human_format(dict_size * sizeof(*A), hdsize);
    
    // Seul le rang 0 affiche l'info pour ne pas spammer la console
    if (rank == 0) {
        printf("Total Dictionary size distributed over %d procs.\n", nprocs);
        printf("Local Dictionary size: %sB per node\n", hdsize);
    }

    A = malloc(sizeof(*A) * dict_size);
    if (A == NULL)
        err(1, "impossible to allocate the dictionnary");
        
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
int golden_claw_search(int maxres, u64 k1[], u64 k2[])
{
    double start = wtime();
    u64 N = 1ull << n;
    for (u64 x = 0; x < N; x++) {
        u64 z = f(x);
        dict_insert(z, x);
    }

    double mid = wtime();
    printf("Fill: %.1fs\n", mid - start);
    
    int nres = 0;
    u64 ncandidates = 0;
    u64 x[256];
    for (u64 z = 0; z < N; z++) {
        u64 y = g(z);
        int nx = dict_probe(y, 256, x);
        assert(nx >= 0);
        ncandidates += nx;
        for (int i = 0; i < nx; i++)
            if (is_good_pair(x[i], z)) {
            	if (nres == maxres)
            		return -1;
            	k1[nres] = x[i];
            	k2[nres] = z;
            	printf("SOLUTION FOUND!\n");
            	nres += 1;
            }
    }
    printf("Probe: %.1fs. %" PRId64 " candidate pairs tested\n", wtime() - mid, ncandidates);
    return nres;
}



int golden_claw_search2(int maxres, u64 k1[], u64 k2[], MPI_Comm comm)
{
    double start_time = wtime();
    u64 N = 1ull << n;

    int rank, nprocs;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &nprocs);

    // 1. Découpage de l'espace de recherche (Work Partitioning)
    u64 range = N / nprocs;
    u64 remainder = N % nprocs;
    
    // Calcul précis pour gérer le reste si N n'est pas divisible par nprocs
    u64 start_idx = rank * range + (rank < remainder ? rank : remainder);
    u64 my_count = range + (rank < remainder ? 1 : 0);
    u64 end_idx = start_idx + my_count;

    // --- PHASE 1 : REMPLISSAGE (FILL) ---
    // On parcourt NOS x, on calcule f(x), mais on ne stocke pas encore
    // car le stockage dépendra du hash (étape suivante)

    // `send_buffers` est un tableau de tableaux. 
    // send_buffers[i] contiendra la liste des données à envoyer au processus 'i'.
    entry_data **send_buffers = malloc(nprocs * sizeof(entry_data*));

    // `send_counts` retient combien d'éléments on a déjà mis dans chaque buffer.
    int *send_counts = calloc(nprocs, sizeof(int));

    // `send_capacities` retient la taille mémoire actuelle de chaque buffer.
    int *send_capacities = malloc(nprocs * sizeof(int));

    // Estimation de la taille : Si le hash est uniforme, chaque processus recevra environ 1/nprocs des données.
    // On prend une marge de sécurité de 20% (* 1.2) pour éviter de faire trop de reallocs au début.
    // On ajoute +128 pour être sûr d'avoir un minimum si my_count est petit.
    int estimated_size = (my_count / nprocs) * 1.2 + 128;

    for(int i = 0; i < nprocs; i++) {
        send_capacities[i] = estimated_size;
        // On alloue la mémoire pour les données
        send_buffers[i] = malloc(send_capacities[i] * sizeof(entry_data));
        
        if (send_buffers[i] == NULL) {
            fprintf(stderr, "[Rank %d] Erreur malloc buffer %d\n", rank, i);
            MPI_Abort(comm, 1);
        }
    }



    for (u64 x = start_idx; x < end_idx; x++) {
        u64 z = f(x);
        
        u64 hash = murmur64(z);
        int target = hash % nprocs; // Le numéro du processus destinataire
        
        // Vérification de la place dans le buffer
        if (send_counts[target] >= send_capacities[target]) {
            // Le bac est plein ! On l'agrandit (x2).
            send_capacities[target] *= 2;
            entry_data *temp = realloc(send_buffers[target], send_capacities[target] * sizeof(entry_data));
            
            if (temp == NULL) {
                fprintf(stderr, "[Rank %d] Erreur realloc buffer %d\n", rank, target);
                MPI_Abort(comm, 1);
            }
            send_buffers[target] = temp;
        }

        // On stocke la donnée dans le buffer
        // Note : On stocke {z, x} car z est la clé de recherche (hash) et x la valeur retrouvée.
        send_buffers[target][send_counts[target]].data_val = z;   // La clé f(x)
        send_buffers[target][send_counts[target]].candidate = x;  // La valeur x
        
        send_counts[target]++; // On incrémente le compteur
        // TODO (Etape 2): 
        // 1. Calculer target = murmur64(z) % nprocs
        // 2. Stocker (z, x) dans un buffer d'envoi pour 'target'
    }

    double mid_time = wtime();
    // Petit print de debug (seulement le rank 0 pour pas spammer)
    if (rank == 0) printf("Fill preparation done in %.2fs. Ready to exchange.\n", mid_time - start_time);

    // TODO (Etape 3): Échanger les buffers (MPI_Alltoall) et insérer dans le dict local

    // 3.1 Échanger les tailles (Combien vais-je recevoir de chacun ?)
    int *recv_counts = malloc(nprocs * sizeof(int));

    // Tout le monde envoie son tableau 'send_counts' à tout le monde.
    // À la fin, recv_counts[i] contiendra le nombre d'éléments que le processus 'i' va m'envoyer.
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    // 3.2 Préparer les déplacements (displacements) pour MPI_Alltoallv
    // MPI a besoin de savoir où commence les données de chaque processus dans le buffer contigu.
    int *send_displs = malloc(nprocs * sizeof(int));
    int *recv_displs = malloc(nprocs * sizeof(int));

    u64 total_send = 0;
    u64 total_recv = 0;
    
    // Calcul des offsets (déplacements)
    send_displs[0] = 0;
    recv_displs[0] = 0;
    total_send = send_counts[0];
    total_recv = recv_counts[0];
    
    for (int i = 1; i < nprocs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        
        total_send += send_counts[i];
        total_recv += recv_counts[i];
    }

    // 3.3 Aplatir les buffers d'envoi (Flattening)
    // MPI_Alltoallv veut un seul gros tableau contigu pour l'envoi.
    // Nous devons copier nos petits buffers dispersés (send_buffers[i]) dans un gros (flat_send_buf).
    
    entry_data *flat_send_buf = malloc(total_send * sizeof(entry_data));
    
    for (int i = 0; i < nprocs; i++) {
        // Copie mémoire rapide
        // Dest : &flat_send_buf[offset]
        // Source : send_buffers[i]
        // Taille : send_counts[i] * sizeof(entry_data)
        if (send_counts[i] > 0) {
             // memcpy est plus rapide qu'une boucle for
             memcpy(&flat_send_buf[send_displs[i]], send_buffers[i], send_counts[i] * sizeof(entry_data));
        }
        // On peut libérer les petits buffers individuels maintenant pour gagner de la RAM
        free(send_buffers[i]); 
    }
    free(send_buffers); // Libération du tableau de pointeurs

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
    free(send_counts); free(recv_counts);
    free(send_displs); free(recv_displs);
    free(send_capacities);
    // On garde les buffers _b pour la phase suivante ou on les free et realloc ? 
    // Mieux vaut tout free pour être propre, on réallouera pour la phase Probe.
    free(s_cnts_b); free(r_cnts_b); free(s_displs_b); free(r_displs_b);

    MPI_Barrier(comm); // Attendre que tout le monde ait fini de remplir
    mid_time = wtime();
    if (rank == 0) printf("Fill phase logic done in %.1fs\n", mid_time - start_time);


    // --- PHASE 2 : SONDAGE (PROBE) ---
    // Réallocation des structures pour la phase 2 (on repart à zéro)
    
    send_buffers = malloc(nprocs * sizeof(entry_data*));
    send_counts = calloc(nprocs, sizeof(int));
    send_capacities = malloc(nprocs * sizeof(int));
    
    // On réutilise la même estimation de taille
    estimated_size = (my_count / nprocs) * 1.2 + 128;

    for(int i = 0; i < nprocs; i++) {
        send_capacities[i] = estimated_size;
        send_buffers[i] = malloc(send_capacities[i] * sizeof(entry_data));
        if (send_buffers[i] == NULL) MPI_Abort(comm, 1);
    }

    // --- 4. Remplissage des buffers avec g(y) ---
    // On parcourt NOS y (les mêmes indices que pour x)
    for (u64 y_val = start_idx; y_val < end_idx; y_val++) {
        u64 z_prime = g(y_val); 

        // Règle d'or : On utilise EXACTEMENT le même hachage que pour f(x)
        // Si f(x) == g(y), alors murmur64(f(x)) == murmur64(g(y))
        // Donc ils iront vers le même 'target'.
        u64 hash = murmur64(z_prime);
        int target = hash % nprocs;

        // Gestion capacité (identique phase 1)
        if (send_counts[target] >= send_capacities[target]) {
            send_capacities[target] *= 2;
            entry_data *temp = realloc(send_buffers[target], send_capacities[target] * sizeof(entry_data));
            if (temp == NULL) MPI_Abort(comm, 1);
            send_buffers[target] = temp;
        }

        // On stocke la requête
        // data_val = z_prime (la valeur qu'on cherche dans le dico)
        // candidate = y_val (la valeur k2 potentielle)
        send_buffers[target][send_counts[target]].data_val = z_prime;
        send_buffers[target][send_counts[target]].candidate = y_val;
        send_counts[target]++;
    }

    // --- 5. Échange MPI et Vérification (Probe) ---

    // 5.1 Échange des tailles
    // On doit réallouer recv_counts car on l'avait free
    recv_counts = malloc(nprocs * sizeof(int));
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm);

    // 5.2 Recalcul des offsets
    send_displs = malloc(nprocs * sizeof(int));
    recv_displs = malloc(nprocs * sizeof(int));
    
    send_displs[0] = 0; recv_displs[0] = 0;
    total_send = send_counts[0];
    total_recv = recv_counts[0];

    for (int i = 1; i < nprocs; i++) {
        send_displs[i] = send_displs[i-1] + send_counts[i-1];
        recv_displs[i] = recv_displs[i-1] + recv_counts[i-1];
        total_send += send_counts[i];
        total_recv += recv_counts[i];
    }

    // 5.3 Aplatir (Flatten)
    flat_send_buf = malloc(total_send * sizeof(entry_data));
    for (int i = 0; i < nprocs; i++) {
        if (send_counts[i] > 0) {
             memcpy(&flat_send_buf[send_displs[i]], send_buffers[i], send_counts[i] * sizeof(entry_data));
        }
        free(send_buffers[i]); 
    }
    free(send_buffers);

    // 5.4 Buffer de réception
    recv_buf = malloc(total_recv * sizeof(entry_data));
    
    // 5.5 Conversion Bytes
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

    // 5.6 L'échange (Probe Data)
    MPI_Alltoallv(flat_send_buf, s_cnts_b, s_displs_b, MPI_BYTE,
                  recv_buf,      r_cnts_b, r_displs_b, MPI_BYTE,
                  comm);

    // 5.7 Vérification Locale (Le cœur de l'attaque)
    
    int nres = 0;       // Solutions trouvées localement
    u64 x_candidates[256]; // Buffer pour récupérer les collisions
    
    for (u64 i = 0; i < total_recv; i++) {
        u64 val_to_find = recv_buf[i].data_val;  // g(y)
        u64 candidate_y = recv_buf[i].candidate; // y

        // On regarde dans NOTRE dictionnaire local si on a ce hash
        int nx = dict_probe(val_to_find, 256, x_candidates);
        
        // dict_probe remplit x_candidates avec les x tels que f(x) == val_to_find
        for (int k = 0; k < nx; k++) {
            // On a une collision hash ! 
            // x_candidates[k] vient du dico (donc c'est un k1 potentiel)
            // candidate_y vient du réseau (donc c'est un k2 potentiel)
            
            if (is_good_pair(x_candidates[k], candidate_y)) {
                if (nres < maxres) {
                    k1[nres] = x_candidates[k];
                    k2[nres] = candidate_y;
                    nres++;
                    printf("[Rank %d] FOUND SOLUTION: k1=%lx, k2=%lx\n", rank, x_candidates[k], candidate_y);
                }
            }
        }
    }

    // --- Nettoyage Final ---
    free(flat_send_buf); free(recv_buf);
    free(send_counts); free(recv_counts);
    free(send_displs); free(recv_displs);
    free(send_capacities);
    free(s_cnts_b); free(r_cnts_b); free(s_displs_b); free(r_displs_b);

    // Note: Chaque processus retourne SON nombre de solutions trouvées.
    // Le main ne verra que les solutions du rank qui a trouvé (via les printf).
    // Pour être 100% propre, il faudrait faire un MPI_Gather des solutions vers le rank 0,
    // mais pour ce projet, le printf suffit généralement.
    
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

    // Seul le rank 0 affiche les infos de démarrage
    if (rank == 0) {
        printf("Running with n=%d, C0=(%08x, %08x) and C1=(%08x, %08x)\n", 
            (int) n, C[0][0], C[0][1], C[1][0], C[1][1]);
    }

    // Allocation distribuée
    dict_setup_mpi(1.125 * (1ull << n), rank, nprocs);

    /* search */
    u64 k1[16], k2[16];
    
    // Appel de la recherche
    int nkey = golden_claw_search2(16, k1, k2, MPI_COMM_WORLD);
    
    // --- CORRECTION ICI ---
    // On enlève assert(nkey > 0) car tous les processus ne trouvent pas forcément une solution.
    
    // On peut faire une barrière pour que l'affichage reste propre (optionnel)
    MPI_Barrier(MPI_COMM_WORLD);

    /* validation */
    // La boucle ne s'exécute QUE si nkey > 0 (donc seulement sur le processus gagnant)
    for (int i = 0; i < nkey; i++) {
        // Double vérification locale
        assert(f(k1[i]) == g(k2[i]));
        assert(is_good_pair(k1[i], k2[i]));     
        
        // On affiche fièrement le résultat
        printf("#################################################################\n");
        printf("[Rank %d] FINAL VALIDATED SOLUTION: (%" PRIx64 ", %" PRIx64 ")\n", rank, k1[i], k2[i]);
        printf("#################################################################\n");
    }

    // Nettoyage final
    free(A);
    MPI_Finalize();
    
    return 0;
}