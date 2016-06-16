#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <mpi.h>

/////////////////////////////////////////////////////////////
//	LIST DATA STRUCTURE	/////////////////////////////
/////////////////////////////////////////////////////////////

typedef struct node_t{
	int data;
	struct node_t *next;
}node;

typedef struct list_t{
	node *Head;
}list;

node *new_node(int data1){
	node *p;
	p = (node *) malloc(sizeof(node));
	p->next = NULL;
	p->data = data1;
	return p;
}

list *new_list(){
	list *list1;
	list1 = (list *) malloc(sizeof(list));
	list1->Head = NULL;
	return list1;
}

list *add_to_head(list *list1, node *node1){
	if (list1->Head == NULL){	
		list1->Head = node1;
	}
	else {
		node1->next = list1->Head;
		list1->Head = node1;
	}
	return list1;
}

////////////////////////////////////////////////////////////
//////////////////////// Block and proc Structure //////////
////////////////////////////////////////////////////////////

typedef struct block_struct_t{
	int block_id;
	int block_n;
	int block_nnz;
	int block_dim;
	int *local_to_global_idx;
	int *global_to_local_idx;
	int *block_rowp;
	int *block_coli;
	double *block_val;
	double *x;
	double *sum;

	int *dsend;
	int *drecv;

	int *sendbuffer;
	int *recvbuffer;

	int bsend;
	int brecv;
	//int psend;
	//int precv;

	bool *flag;
	int *perm;
	//int *reindep;
	bool *copyflag;
	int block_counter;
	int no_ind;
	int *idx_2_order;
	int *last_xi;
	list **list_dep;
	list **list_ind;
	int *num_dep_diag;
	int *num_dep_offdiag;

	bool *dual_check;

}block_struct;

typedef struct proc_struct_t{
	int nbpp;
	block_struct **block_list;
}proc_struct;

block_struct *new_block(){
	block_struct *block1;
	block1 = (block_struct *) calloc(1,sizeof(block_struct));
	block1->block_id = -1;
	block1->block_n = 0;
	block1->block_nnz = 0;
	block1->block_dim = 0;
	block1->local_to_global_idx = NULL;
	block1->global_to_local_idx = NULL;
	block1->block_rowp = NULL;
	block1->block_coli = NULL;
	block1->block_val = NULL;
	block1->x = NULL;
	block1->sum = NULL;

	block1->dsend = NULL;
	block1->drecv = NULL;
	block1->sendbuffer = NULL;
	block1->recvbuffer = NULL;
	block1->bsend = MPI_PROC_NULL;
	block1->brecv = MPI_PROC_NULL;
	//block1->psend = MPI_PROC_NULL;
	//block1->precv = MPI_PROC_NULL;

	block1->flag = NULL;
	block1->perm = NULL;
	//block1->reindep = NULL;
	block1->copyflag = NULL;
	
	block1->block_counter = 0;
	block1->no_ind = 0;

	block1->last_xi = NULL;
	block1->idx_2_order = NULL;
	block1->list_dep = NULL;
	block1->list_ind = NULL;
	
	block1->num_dep_diag = NULL;
	block1->num_dep_offdiag = NULL;

	return block1;
}

proc_struct *new_proc(int nbpp1){
	proc_struct *proc1;
	proc1 = (proc_struct *) calloc(1,sizeof(proc_struct));
	proc1->nbpp = nbpp1;
	proc1->block_list =(block_struct **) calloc(nbpp1,sizeof(block_struct *));
	return proc1;
}

///////////////////////////////////////////////////////////////////////////////////////////


void process_row_ind(proc_struct *proc1, int r){
	int k1,k2,k3,k4,k5,k6,k7,i1,i2,i3,i4;
	double buf;
	node *p;
	MPI_Request Req1;	
	
	// Solving for X[i]
	k1 = proc1->block_list[0]->block_rowp[r];
	k2 = proc1->block_list[0]->block_rowp[r+1]-1;
	for (i1=k1; i1<k2; i1++){
		i2 = proc1->block_list[0]->block_coli[i1];
		proc1->block_list[0]->sum[r] += proc1->block_list[0]->x[i2];
	}
	proc1->block_list[0]->x[r] = 1 - proc1->block_list[0]->sum[r];
	//printf ("proc:%d:: x[%d]: %lf\n", proc1->block_list[0]->block_id ,r, proc1->block_list[0]->x[r]);

	// Checking the list for off diagonal rows
	k3 = proc1->block_list[0]->idx_2_order[r];
	p = proc1->block_list[0]->list_ind[k3]->Head;
	while (p != NULL){
		i1 = p->data; // Actual local row iindex
		k4 = proc1->block_list[0]->block_rowp[i1];
		k5 = proc1->block_list[0]->block_rowp[i1+1];
		for (i2=k4; i2<k5; i2++){
			i3 = proc1->block_list[0]->block_coli[i2];
			proc1->block_list[0]->sum[i1] += proc1->block_list[0]->x[i3];
		}
		proc1->block_list[0]->num_dep_offdiag[i1-proc1->block_list[0]->block_dim] = -1;
		if (proc1->block_list[0]->dual_check[i1] == true){
			k6 = proc1->block_list[0]->local_to_global_idx[i1];	// tag
			k1 = proc1->block_list[0]->block_dim;
			k7 = proc1->block_list[0]->dsend[i1-k1];	// destination
			buf = proc1->block_list[0]->sum[i1];	// buffer
			//MPI_Isend(&buf, 1, MPI_DOUBLE, k7, k6, MPI_COMM_WORLD, &Req1); // Check Isend
			MPI_Send(&buf, 1, MPI_DOUBLE, k7, k6, MPI_COMM_WORLD);
			printf (" %d is sending to %d, for %d row.\n", proc1->block_list[0]->block_id, k7, k6);
			
		}
		p = p->next;
	}			
}


void process_row_dep(proc_struct *proc1, int r){
	int k1,k2,k3,k4,k5,k6,k7,i1,i2,i3,i4, i5;
	double buf;
	node *p;
	MPI_Request Req1;	
	
	// Solving for X[i]
	k1 = proc1->block_list[0]->block_rowp[r];
	k2 = proc1->block_list[0]->block_rowp[r+1]-1;
	for (i1=k1; i1<k2; i1++){
		i2 = proc1->block_list[0]->block_coli[i1];
		proc1->block_list[0]->sum[r] += proc1->block_list[0]->x[i2];
	}
	proc1->block_list[0]->x[r] = 1 - proc1->block_list[0]->sum[r];
	//printf ("proc %d:: x[%d]: %lf\n",proc1->block_list[0]->block_id, r, proc1->block_list[0]->x[r]);

	// Checking the dependent row on this row
	k3 = proc1->block_list[0]->idx_2_order[r] - proc1->block_list[0]->no_ind;
	p = proc1->block_list[0]->list_dep[k3]->Head;
	while (p != NULL){
		i1 = p->data; // Actual local row iindex
		i2 = proc1->block_list[0]->idx_2_order[i1];	// order
 
		if  (i1 >= proc1->block_list[0]->block_dim){	// num_dep_offdiag
			i3 = i1 - proc1->block_list[0]->block_dim;
			proc1->block_list[0]->num_dep_offdiag[i3] = proc1->block_list[0]->num_dep_offdiag[i3] - 1;

			if (proc1->block_list[0]->num_dep_offdiag[i3] == 0){
				k4 = proc1->block_list[0]->block_rowp[i1];
				k5 = proc1->block_list[0]->block_rowp[i1+1];
				for (i4=k4; i4<k5; i4++){
					i5 = proc1->block_list[0]->block_coli[i4];
					proc1->block_list[0]->sum[i1] += proc1->block_list[0]->x[i5];
				}
				proc1->block_list[0]->num_dep_offdiag[i3] = proc1->block_list[0]->num_dep_offdiag[i3] - 1;
				if (proc1->block_list[0]->dual_check[i1] == true){	// send

					//proc1->block_list[0]->x[i1] = 1 - proc1->block_list[0]->sum[i1];
					k6 = proc1->block_list[0]->local_to_global_idx[i1];	// tag
					k1 = proc1->block_list[0]->block_dim;
					k7 = proc1->block_list[0]->dsend[i1-k1];	// destination
					buf = proc1->block_list[0]->sum[i1];
					//MPI_Isend(&buf, 1, MPI_DOUBLE, k7, k6, MPI_COMM_WORLD, &Req1);
					MPI_Send(&buf, 1, MPI_DOUBLE, k7, k6, MPI_COMM_WORLD);
					printf ("%d is sending to %d, for %d row..\n",proc1->block_list[0]->block_id,k7,k6);
					proc1->block_list[0]->num_dep_offdiag[i3] = proc1->block_list[0]->num_dep_offdiag[i3] - 1;
				}
			}						
			//p = p->next;
		}

		else if (i2 >= proc1->block_list[0]->no_ind){	// Num_dep_diag
			i3 = i2 - proc1->block_list[0]->no_ind;	// among num_dep_diag
			proc1->block_list[0]->num_dep_diag[i3] = proc1->block_list[0]->num_dep_diag[i3] - 1;
			//printf ("Main hun %d  = %d .... /n", proc1->block_list[0]->local_to_global_idx[i1], proc1->block_list[0]->num_dep_diag[i3]);

			if ((proc1->block_list[0]->num_dep_diag[i3] == 0) && (proc1->block_list[0]->dual_check[i1] == true)){
				//printf ("Main hun %d  of %d.... /n", proc1->block_list[0]->local_to_global_idx[i1], proc1->block_list[0]->block_id);
				process_row_dep(proc1,i1);
			}
		}
	
		else {
			printf ("Error in list for dep");
		}		
		p = p->next;
	}
		
}

void triangular_solve(proc_struct *proc1){
	int k1,k2,k3,k4,k5,k6,k7,i1,i2,i3,i4,j1,j2,j3;
	int myid = proc1->block_list[0]->block_id;
	MPI_Status *status = (MPI_Status *) malloc(sizeof(MPI_Status));
	MPI_Status *status1 = (MPI_Status *) malloc(sizeof(MPI_Status));
	int recv_count = 0;
	double tempsum;

	k2 = proc1->block_list[0]->block_dim;
	k1 = proc1->block_list[0]->block_n;
	k3 = proc1->block_list[0]->no_ind;

	proc1->block_list[0]->x = (double *) calloc(k2, sizeof(double));
	proc1->block_list[0]->sum = (double *) calloc(k1, sizeof(double));
	proc1->block_list[0]->dual_check = (bool *) malloc(k1*sizeof(bool));
	
	if (proc1->block_list[0]->block_id == 0){
		for (i1=0; i1<k1; i1++)	proc1->block_list[0]->dual_check[i1] = true;
	}
	else {
	for (i1=0; i1<k1; i1++){
		if (proc1->block_list[0]->drecv[i1] == -1){
			proc1->block_list[0]->dual_check[i1] = true;
		}
		else{
			 proc1->block_list[0]->dual_check[i1] = false;
			 recv_count += 1;
		}
	}}
	
	//printf ("blockid :%d \n", myid);
	
	
	/// Processing Independent Rows
	for (i1=0; i1<k3; i1++){
		k4 = proc1->block_list[0]->perm[i1]; // local row index in block - diag
		process_row_ind(proc1, k4);
	}
	
	// Check for data to be recieved
	//printf ("%d recv count : %d\n",proc1->block_list[0]->block_id, recv_count);
	/*if (proc1->block_list[0]->block_id == 2){	
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
		printf("proc: %d   :: Source tag: %d %d \n",proc1->block_list[0]->block_id,(*status).MPI_SOURCE, (*status).MPI_TAG);
		MPI_Recv(&tempsum, 1, MPI_DOUBLE, (*status).MPI_SOURCE, (*status).MPI_TAG, MPI_COMM_WORLD, status);
		//MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status);
		//printf("proc: %d   :: Source tag: %d %d \n",proc1->block_list[0]->block_id,(*status).MPI_SOURCE, (*status).MPI_TAG);
		//MPI_Recv(&tempsum, 1, MPI_DOUBLE, (*status).MPI_SOURCE, (*status).MPI_TAG, MPI_COMM_WORLD, status);
	}*/
	//if (proc1->block_list[0]->block_id != 2){
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	for (i1=0; i1<recv_count; i1++){
		//printf("CHUTIYA\n");
		MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status); //
		//printf("proc: %d   :: Source tag: %d %d \n",proc1->block_list[0]->block_id,(*status).MPI_SOURCE, (*status).MPI_TAG);
		i2 = (*status).MPI_TAG;	// Global row index
		i3 = proc1->block_list[0]->global_to_local_idx[i2]; // Local row index
		i4 = (*status).MPI_SOURCE;
		//printf("proc: %d   :: Source tag: %d %d \n",proc1->block_list[0]->block_id,(*status).MPI_SOURCE, i2);
		MPI_Recv(&tempsum, 1, MPI_DOUBLE, i4, i2, MPI_COMM_WORLD, status1);
		proc1->block_list[0]->sum[i3] = proc1->block_list[0]->sum[i3] + tempsum;

		
		// Recieved massage for dep-diag
		if (i3 < proc1->block_list[0]->block_dim){
			proc1->block_list[0]->dual_check[i3] = true;
			j1 = proc1->block_list[0]->idx_2_order[i3] - proc1->block_list[0]->no_ind;
			j2 = proc1->block_list[0]->num_dep_diag[j1];
			//printf ("***Dependent proc %d:: %d: %d\n",proc1->block_list[0]->block_id, i3, j2);
			if (j2 == 0){
					process_row_dep(proc1, i3);
					//proc1->block_list[0]->num_dep_diag[i] = proc1->block_list[0]->num_dep_diag[i3] - 1;
					//printf ("***Dependent proc %d:: x[%d]: %lf\n",proc1->block_list[0]->block_id, i3, proc1->block_list[0]->x[i3]);		
			}	
		}

		// Recieved massage for off-diag
		else{
			if (proc1->block_list[0]->num_dep_offdiag[i3-proc1->block_list[0]->block_dim] == -1){
				k7 = proc1->block_list[0]->dsend[i3-proc1->block_list[0]->block_dim];	// destination
				MPI_Send(&proc1->block_list[0]->sum[i3], 1, MPI_DOUBLE, k7, i2, MPI_COMM_WORLD);
			}
			else{
				proc1->block_list[0]->dual_check[i3] = true;
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}


void find_last_xi(proc_struct *proc1){
	int k1,k2,k3,k4,k5,k6,k7,i1,i2,i3,i4;
	k1 = proc1->block_list[0]->block_dim;
	k2 = proc1->block_list[0]->block_n;

	proc1->block_list[0]->idx_2_order = (int *) calloc(k1, sizeof(int));
	for (i1=0; i1<k1; i1++)	proc1->block_list[0]->idx_2_order[proc1->block_list[0]->perm[i1]] = i1;

	proc1->block_list[0]->last_xi = (int *) calloc(k2-k1, sizeof(int));
	for (i1=0; i1<(k2-k1); i1++)	proc1->block_list[0]->last_xi[i1] = -1;	
	
	for (i1=k1; i1<k2; i1++){
		k3 = proc1->block_list[0]->block_rowp[i1];
		k4 = proc1->block_list[0]->block_rowp[i1+1];
		k5 = -1;
		for (i2=k3; i2<k4; i2++){
			k6 = proc1->block_list[0]->block_coli[i2];
			k7 = proc1->block_list[0]->idx_2_order[k6];
			if (k7 > k5){
				k5 = k7;
				proc1->block_list[0]->last_xi[i1-k1] = k6;
			}	
		}	
	}
	
	/*if (proc1->block_list[0]->block_id == 0){
		for (i1=0; i1<(k2-k1); i1++)	printf("%d ",proc1->block_list[0]->last_xi[i1]);
	}
	if (proc1->block_list[0]->block_id == 0){
		printf("\ndsend next\n");
		for (i1=0; i1<k2-k1; i1++)	printf("%d ",proc1->block_list[0]->dsend[i1]);
	}*/
}

void analys_offdiag(proc_struct *proc1){
	int k1,k2,k3,k4,k5,k6,k7,k8,i1,i2,i3,i4;
	int maxidx = -10, break1, break2;
	node *node1;

	k1 = proc1->block_list[0]->block_dim;
	k2 = proc1->block_list[0]->block_n;
	k3 = proc1->block_list[0]->no_ind;

	if (k3 != 0){
	proc1->block_list[0]->list_ind = (list **) malloc(k3*sizeof(list *));
	for (i1=0; i1<k3; i1++)	proc1->block_list[0]->list_ind[i1] = new_list();
	}
	proc1->block_list[0]->num_dep_offdiag = (int *) calloc(k2-k1, sizeof(int));
		
	for (i1=k1; i1<k2; i1++){
		k4 = proc1->block_list[0]->block_rowp[i1];
		k5 = proc1->block_list[0]->block_rowp[i1+1];
		k6 = -1;
		break2 = -5;
		break1 = -5; // No break
		for (i2=k4; i2<k5; i2++){
			break2 = +5;
			k7 = proc1->block_list[0]->block_coli[i2];
			k8 = proc1->block_list[0]->idx_2_order[k7];
			if (k8 >= k3){
				break1 = i2;	
				break;
			}
			if (k8 > k6){
				k6 = k8;
				maxidx = k8;
			}			
		}
		//printf ("%d : break1 %d, break2 %d, maxidx %d\n", proc1->block_list[0]->block_id, break1, break2, maxidx);
		if ((break1 == -5) && (break2 == +5)){
			node1 = new_node(i1);
			proc1->block_list[0]->list_ind[maxidx] = add_to_head(proc1->block_list[0]->list_ind[maxidx], node1);
		}
		if ((break1 != -5) && (break2 == +5)){
			for (i2=break1; i2<k5; i2++){
				k7 = proc1->block_list[0]->block_coli[i2];
				k8 = proc1->block_list[0]->idx_2_order[k7];
				if (k8 >= k3){
					proc1->block_list[0]->num_dep_offdiag[i1-k1] += 1;
					node1 = new_node(i1);
					proc1->block_list[0]->list_dep[k8-k3] = add_to_head(proc1->block_list[0]->list_dep[k8-k3], node1);
				}
			}
		}
	}	
}


void analys_dep_diag(proc_struct *proc1){
	int k1,k2,k3,k4,k5,k6,k7,i1,i2,i3,i4;
	node *node1;

	k1 = proc1->block_list[0]->block_dim;
	k2 = proc1->block_list[0]->block_id;
	k3 = proc1->block_list[0]->no_ind;
	//printf ("%d : %d\n", k2, k3);

	if (k1 != k3){
	proc1->block_list[0]->num_dep_diag = (int *) calloc(k1-k3, sizeof(int));
	proc1->block_list[0]->list_dep = (list **) malloc((k1-k3)*sizeof(list *));
	for (i1=0; i1<(k1-k3); i1++)	proc1->block_list[0]->list_dep[i1] = new_list();

	//printf("Malloc done\n");

	for (i1=k3; i1<k1; i1++){
		i2 = proc1->block_list[0]->perm[i1];	// Actual local row index
		//printf ("%d, %d : %d \n", k2, i1, i2);
		k4 = proc1->block_list[0]->block_rowp[i2];
		k5 = proc1->block_list[0]->block_rowp[i2+1]-1;
		for (i3=k4; i3<k5; i3++){
			k6 = proc1->block_list[0]->block_coli[i3];
			//printf ("%d k6 %d\n",k2, proc1->block_list[0]->idx_2_order[0]);
			k7 = (proc1->block_list[0]->idx_2_order[k6]) - k3;
			
			if (proc1->block_list[0]->flag[k6] == true){
				proc1->block_list[0]->num_dep_diag[i1-k3] += 1;
				node1 = new_node(i2);
				proc1->block_list[0]->list_dep[k7] = add_to_head(proc1->block_list[0]->list_dep[k7], node1);
			}
		}
			
	}
	}
	/*if (k2 == 1){
		node *p = proc1->block_list[0]->list_dep[0]->Head;
		while (p != NULL){	
			printf("%d ",p->data);
			p = p->next;
		}
	}*/
}


void copyRowRecursive(proc_struct *proc1, int r){
	int k1,k2,k3,k4,k5,i1,i2,i3,i4;

	k1 = proc1->block_list[0]->block_dim;
	k2 = proc1->block_list[0]->block_id;	
	k3 = proc1->block_list[0]->block_rowp[r];
        k4 = proc1->block_list[0]->block_rowp[r+1]-2;

	for (i1=k4; i1>=k3; i1--){
		k5 = proc1->block_list[0]->block_coli[i1];
		//printf("%d %d:%d %d\n",k2, r, k5, proc1->block_list[0]->block_counter);
		if (proc1->block_list[0]->copyflag[k5] == false)	copyRowRecursive(proc1, k5);
	}
	
	proc1->block_list[0]->perm[proc1->block_list[0]->block_counter] = r;
        proc1->block_list[0]->copyflag[r] = true;
        proc1->block_list[0]->block_counter += 1;
	//printf("%d %d: %d\n",k2, r, proc1->block_list[0]->block_counter);
}


void reorder_dig (proc_struct *proc1){
	int i1,i2,i3, k1, k2, k3, k4, k5;

	k1 = proc1->block_list[0]->block_dim;
	k2 = proc1->block_list[0]->block_id;

	proc1->block_list[0]->flag = (bool *) malloc(sizeof(bool)*k1);

	for (i1=0; i1<k1; i1++)	proc1->block_list[0]->flag[i1] = false; // Default all independent

	if (k2 != 0){
		for (i1=0; i1<k1; i1++){ // For each row in diag
			if (proc1->block_list[0]->drecv[i1] != -1){
				proc1->block_list[0]->flag[i1] = true; // dependent
			}
			else {
				k3 = proc1->block_list[0]->block_rowp[i1];
				k4 = proc1->block_list[0]->block_rowp[i1+1];
				for (i2=k3; i2<k4; i2++){
					k5 = proc1->block_list[0]->block_coli[i2];
					if (proc1->block_list[0]->flag[k5] == true){
						proc1->block_list[0]->flag[i1] = true;
						break;
					}
				}
			}
		}
	}

//	printf("COLI Check ...........\n");
//	for (i1=0; i1<proc1->block_list[0]->block_nnz; i1++)	printf("%d:%d \n",k2, proc1->block_list[0]->block_coli[i1]);
		
//	}

//	for (i1=0; i1<k1; i1++)	printf("%d: %d\n", k2, proc1->block_list[0]->flag[i1]);

////////////////
	proc1->block_list[0]->copyflag = (bool *) malloc(sizeof(bool)*k1);
	proc1->block_list[0]->perm = (int *) calloc(k1,sizeof(int));

	for (i1=0; i1<k1; i1++)	proc1->block_list[0]->copyflag[i1] = false; // not copied

	for (i1=k1-1; i1>=0; i1--){	// Rows in backward direction
		if (proc1->block_list[0]->flag[i1] == false && proc1->block_list[0]->copyflag[i1] == false){
			copyRowRecursive(proc1, i1);
		}
	}

	proc1->block_list[0]->no_ind = proc1->block_list[0]->block_counter;
	for (i1=0; i1<k1; i1++){     // Rows in forward direction
		if (proc1->block_list[0]->flag[i1] == true){
			proc1->block_list[0]->perm[proc1->block_list[0]->block_counter] = i1;
			proc1->block_list[0]->copyflag[i1] = true;
			proc1->block_list[0]->block_counter += 1;
		}
	}
	//if(k2 == 0){
	//printf("\n");
	//for (i1=0; i1<k1; i1++)	printf("%d ",proc1->block_list[0]->perm[i1]);
	//}
///////////////////////////////////////////////////////////	
	proc1->block_list[0]->idx_2_order = (int *) calloc(k1, sizeof(int));
	for (i1=0; i1<k1; i1++)	proc1->block_list[0]->idx_2_order[proc1->block_list[0]->perm[i1]] = i1;
}

int main(int argc, char **argv)
{

	// Variables
	int nproc, proc_id;
	int Gblock_dim= atoi(argv[1]);
	int  nblock;	// nblock = no. of blocks
	int mat_dim=9, mat_nnz=14 , Gnbpp;	// nbpp = no. of blocks per process

	MPI_Status stat;

	int i1,i2,i3,i4,j1,j2,j3,k1,k2,k3;

	// MPI Initialization
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	nblock = mat_dim/Gblock_dim;
	Gnbpp = nblock/nproc;

	proc_struct *proc1;
	proc1 = new_proc(Gnbpp);
	block_struct *block1;

	int myblock_id[Gnbpp], myn[Gnbpp], mynnz[Gnbpp];
	int endIdx[Gnbpp], startIdx[Gnbpp], counter[Gnbpp];
	//int *block_rowp[nbpp], *block_coli[nbpp];

	for (i1=0; i1<Gnbpp; i1++)	myblock_id[i1] = (i1*nproc) + proc_id;
	for (i1=0; i1<Gnbpp; i1++){
		startIdx[i1] = myblock_id[i1]*Gblock_dim;
		endIdx[i1] = startIdx[i1] + Gblock_dim - 1;
		mynnz[i1] = 0;
		counter[i1] = 0;
	}

/*
	READING MATRIX DATA IN LOCAL CSR FORMATE FOR EACH BLOCK, GIVEN CSC Global
*/
	FILE *fp1, *fp2, *fp3;
	char skip1[100];
	
	for (i1=0; i1<Gnbpp; i1++){
		myn[i1] = mat_dim - (myblock_id[i1]*Gblock_dim);
		//block_rowp[i1] = (int *) malloc(myn[i1]+1);
		//for (i2=0; i2<=myn[i1]; i2++)	block_rowp[i1][i2] = 0;
	}
	
	fp2 = fopen("coli.dat","r");
	for (i1=0; i1<mat_nnz; i1++){
		fscanf(fp2, "%d\n",&k1);
		for (i2=0; i2<Gnbpp; i2++){
			if (k1 <= endIdx[i2] && k1>=startIdx[i2]){
				mynnz[i2] += 1;
				break;
			}
		}
	}
	
	for (i1=0; i1<Gnbpp; i1++){	//block_coli[i1] = (int *) malloc(mynnz[i1]);
		//printf("%d %d: %d %d, %d %d\n",proc_id,myblock_id[i1],myn[i1],mynnz[i1], startIdx[i1], endIdx[i1]);
	}


	for (i1=0; i1<Gnbpp; i1++){
                block1 = new_block();
		block1->block_id = myblock_id[i1];
                block1->block_n = myn[i1];
		block1->block_nnz = mynnz[i1];
		block1->block_dim = Gblock_dim;
		
		block1->local_to_global_idx = (int *) calloc(myn[i1], sizeof(int));
		for (i2=0; i2<myn[i1]; i2++)	block1->local_to_global_idx[i2] = i2 + startIdx[i1];

		block1->global_to_local_idx = (int *) calloc(myn[i1], sizeof(int));
		for (i2=0; i2<myn[i1]; i2++)	block1->global_to_local_idx[i2] = i2 - startIdx[i1];

		block1->dsend = (int *) calloc((myn[i1]-Gblock_dim), sizeof(int));
		block1->drecv = (int *) calloc(myn[i1], sizeof(int));

		block1->recvbuffer = (int *) calloc(2*(myn[i1]-Gblock_dim), sizeof(int));
		block1->sendbuffer = (int *) calloc(2*myn[i1], sizeof(int));	
		
		if ((myblock_id[i1] != 0) && (myblock_id[i1] != nblock-1)){
			block1->bsend = block1->block_id + 1;
			block1->brecv = block1->block_id - 1;
		}
		if (myblock_id[i1] == 0){
			block1->bsend = block1->block_id + 1;
                        block1->brecv = MPI_PROC_NULL;
		}
		if (myblock_id[i1] == nblock-1){
			block1->brecv = block1->block_id - 1;
                        block1->bsend = MPI_PROC_NULL;
		}


		block1->block_rowp = (int *) calloc((myn[i1]+1), sizeof(int));
		block1->block_coli = (int *) calloc(mynnz[i1], sizeof(int));
		block1->block_val = (double *) calloc(mynnz[i1], sizeof(double));

		//block1->flag = (bool *) calloc((myn[i1]+1), sizeof(int));
		//block1->block_rowp = (int *) calloc((myn[i1]+1), sizeof(int));
		//block1->block_rowp = (int *) calloc((myn[i1]+1), sizeof(int));

                proc1->block_list[i1] = block1;
        }

	fclose(fp2);
	fp2 = fopen("coli.dat","r");
	fp1 = fopen("rowp.dat","r");
	k1 = 0;

	fgets(skip1, 100, fp1);

	for (i1=0; i1<mat_dim; i1++){
		fscanf(fp1, "%d\n",&k2);
		for (i2=0; i2<(k2-k1); i2++){
			fscanf(fp2, "%d\n",&k3);
			for (i3=0; i3<Gnbpp; i3++){
				if (k3 <= endIdx[i3] && k3>=startIdx[i3]){
					proc1->block_list[i3]->block_rowp[i1+1-startIdx[i3]]  += 1;
					proc1->block_list[i3]->block_coli[counter[i3]++] = k3-startIdx[i3];
					break;
				}
			}	
		}
		k1=k2;
	}

	//MPI_Barrier(MPI_COMM_WORLD);
	for (i1=0; i1<Gnbpp; i1++)
		for (i2=2; i2<=myn[i1]; i2++){
			proc1->block_list[i1]->block_rowp[i2] += proc1->block_list[i1]->block_rowp[i2-1];
		//	printf("%d: %d\n ", myblock_id[i1], proc1->block_list[i1]->block_rowp[i2]);
		}

	
/***************************************************************************************************/
	/*  Dependancy Information/ Data massage */

	for (i1=0; i1<Gnbpp; i1++)
                for (i2=Gblock_dim; i2<myn[i1]; i2++){
                        if (proc1->block_list[i1]->block_rowp[i2] != proc1->block_list[i1]->block_rowp[i2+1])
				proc1->block_list[i1]->dsend[i2-Gblock_dim] = proc_id;
			else  proc1->block_list[i1]->dsend[i2-Gblock_dim] = -1;
                }

	if (proc_id == 0){
		//printf ("send start - %d to %d\n",proc_id,proc1->block_list[0]->bsend);	
		MPI_Send(proc1->block_list[0]->dsend, (myn[0]-Gblock_dim), MPI_INT, proc1->block_list[0]->bsend, 1, MPI_COMM_WORLD);
		//printf ("send end - %d to %d\n",proc_id,proc1->block_list[0]->bsend);


		//printf ("Revs recv start - %d to %d\n",proc_id,proc1->block_list[0]->bsend);
                MPI_Recv(proc1->block_list[0]->recvbuffer, 2*(myn[0]-Gblock_dim), MPI_INT, proc1->block_list[0]->bsend, 2, MPI_COMM_WORLD, &stat);
                //printf ("reverse end - %d to %d\n",proc_id,proc1->block_list[0]->bsend);

		for (i1=0; i1<(myn[0]-Gblock_dim); i1++){
			if (proc1->block_list[0]->recvbuffer[i1] == proc_id){
				proc1->block_list[0]->dsend[i1]=proc1->block_list[0]->recvbuffer[(myn[0]-Gblock_dim)+i1];
			}
			else{
				proc1->block_list[0]->dsend[i1]=-1;
			}	
		}

	}

	else if (proc_id == nproc-1){
		//printf ("recv start - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
		MPI_Recv(proc1->block_list[0]->drecv, myn[0], MPI_INT, proc1->block_list[0]->brecv, 1, MPI_COMM_WORLD, &stat);
		//printf ("recv end - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
		for (i1=0; i1<myn[0]; i1++){
			proc1->block_list[0]->sendbuffer[i1] = proc1->block_list[0]->drecv[i1];
			proc1->block_list[0]->sendbuffer[myn[0]+i1] = proc_id;
		}

		//printf ("reverse send start - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
                MPI_Send(proc1->block_list[0]->sendbuffer, 2*myn[0], MPI_INT, proc1->block_list[0]->brecv, 2, MPI_COMM_WORLD);
                //printf ("send end - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
	}

	else {
		//printf ("send start - %d to %d\n",proc_id,proc1->block_list[0]->bsend);
		MPI_Recv(proc1->block_list[0]->drecv, myn[0], MPI_INT, proc1->block_list[0]->brecv, 1, MPI_COMM_WORLD, &stat);
		//printf ("send end - %d to %d\n",proc_id,proc1->block_list[0]->bsend);
		for (i1=Gblock_dim; i1<myn[0]; i1++){
			k1 = proc1->block_list[0]->drecv[i1];
			k2 = proc1->block_list[0]->dsend[i1-Gblock_dim];
			if ((k1 != -1) && (k2 == -1)){	
				proc1->block_list[0]->dsend[i1-Gblock_dim] = k1;
				proc1->block_list[0]->drecv[i1] = -1;
			}
		}
		//printf ("recv start - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
		MPI_Send(proc1->block_list[0]->dsend, (myn[0]-Gblock_dim), MPI_INT, proc1->block_list[0]->bsend, 1, MPI_COMM_WORLD);
		//printf ("recv end - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
		

		//printf ("Revs recv start - %d to %d\n",proc_id,proc1->block_list[0]->bsend);
                MPI_Recv(proc1->block_list[0]->recvbuffer, 2*(myn[0]-Gblock_dim), MPI_INT, proc1->block_list[0]->bsend, 2, MPI_COMM_WORLD, &stat);
                //printf ("reverse recv end - %d to %d\n",proc_id,proc1->block_list[0]->bsend);
		

		for (i1=0; i1<myn[0]; i1++){
			if (proc1->block_list[0]->drecv[i1] != -1){
				proc1->block_list[0]->sendbuffer[i1+myn[0]] = proc_id;
				proc1->block_list[0]->sendbuffer[i1] = proc1->block_list[0]->drecv[i1];
			}
			else {
				proc1->block_list[0]->sendbuffer[i1+myn[0]] = -1;
				proc1->block_list[0]->sendbuffer[i1] = -1;
			}
		}

		for (i1=0; i1<(myn[0]-Gblock_dim); i1++){
			if (proc1->block_list[0]->recvbuffer[i1] == proc_id){
				proc1->block_list[0]->dsend[i1]=proc1->block_list[0]->recvbuffer[(myn[0]-Gblock_dim)+i1];
			}
			else if (proc1->block_list[0]->recvbuffer[i1] == -1){
				proc1->block_list[0]->dsend[i1]=-1;
			}
			else{
				proc1->block_list[0]->dsend[i1]=-1;
				proc1->block_list[0]->sendbuffer[i1+Gblock_dim] = proc1->block_list[0]->recvbuffer[i1];
				proc1->block_list[0]->sendbuffer[i1+Gblock_dim+myn[0]] = proc1->block_list[0]->recvbuffer[(myn[0]-Gblock_dim)+i1];
			} 	
		}

		//printf ("reverse send start - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
                MPI_Send(proc1->block_list[0]->sendbuffer, 2*myn[0], MPI_INT, proc1->block_list[0]->brecv, 2, MPI_COMM_WORLD);
                //printf ("send end - %d to %d\n",proc_id,proc1->block_list[0]->brecv);
	}
	
	/*for (i2=0; i2<nproc; i2++){
	if (proc_id == 2){
		printf("Process : %d\n", proc_id);
		for (i1=0; i1<myn[0]; i1++)	printf("%d ", proc1->block_list[0]->drecv[i1]);
		printf("\n");
		for (i1=0; i1<(myn[0]-Gblock_dim); i1++) printf("%d ", proc1->block_list[0]->dsend[i1]);
		printf("\n");
	}
	}*/
	
	free (proc1->block_list[0]->sendbuffer);
	free (proc1->block_list[0]->recvbuffer);

	reorder_dig (proc1);
	
	/*if (proc_id == 2){
		printf("Process : %d\n", proc_id);
		for (i1=0; i1<3; i1++)	printf("%d ", proc1->block_list[0]->perm[i1]);
		printf("\n");
		//for (i1=0; i1<3; i1++) printf("%d ", proc1->block_list[0]->flag[i1]);
		//printf("\n");
	}*/

	//find_last_xi(proc1);
	analys_dep_diag(proc1);
	/*if (proc_id == 2){
		printf("Process : %d\n", proc_id);
		for (i1=0; i1<3; i1++)	printf("%d ", proc1->block_list[0]->perm[i1]);
		printf("\n");
		//for (i1=0; i1<3; i1++) printf("%d ", proc1->block_list[0]->flag[i1]);
		//printf("\n");
	}*/

	analys_offdiag(proc1);
	//if (proc_id != 2){
	triangular_solve(proc1);
	//}
	for (i1=0; i1<proc1->block_list[0]->block_dim; i1++){
		printf ("x[%d]: %lf\n",proc1->block_list[0]->local_to_global_idx[i1], proc1->block_list[0]->x[i1]);
	}
/*************************************************************************************************/
// Assuming 1 block per process, block_id = proc_id	
	
//	if (proc_id == 0){
//		for (i1=block)
//		proc1->block_list[0]->
//	}
/*************************************************************************************************/
/* 
	Finalizing the results 
*/

	//printf("Abra ka dabra");
	fclose(fp1);
	fclose(fp2);
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	//printf("YAAFsfjgshfr");
	return 0;
}
