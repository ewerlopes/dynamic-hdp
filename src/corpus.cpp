// part of this code in this file is taken from
// the LDA-C code from prof. David Blei.
// reading the file with lda-c format
//
#include <stdlib.h>
#include <stdio.h>
#include "corpus.h"
#include "utils.h"

void document::free_document()
{
    delete [] words;
    delete [] counts;
    words = NULL;
    counts = NULL;
    length = 0;
    total = 0;
    id = -1;
}

void document::read_data(FILE * fileptr, bool& eof, int _id)
{
    free_document();
    
	int OFFSET = 0;
	int word;
	int count;

	if ((fscanf(fileptr, "%10d", &length) != EOF))
    {
		words = new int [length];
		counts = new int [length];
        for (int n = 0; n < length; n++)
        {
            fscanf(fileptr, "%10d:%10d", &word, &count);
            word = word - OFFSET;
            words[n] = word;
            counts[n] = count;
            total += count;
        }
        id = _id;
		eof = false;
    } else {
		eof = true;
	}
}

int document::get_max_word_id() const
{
	int * index = new int;
	int max_word_id = max(words, length, index);
    delete index;
    return max_word_id;
}



corpus::corpus()
{
    size_vocab = 0;
    total_words = 0;
    num_docs = 0;
}

corpus::~corpus()
{
    for (int i = 0; i < num_docs; i++)
    {
        document * doc = docs[i];
        delete doc;
    }
    docs.clear();

    size_vocab = 0;
    num_docs = 0;
    total_words = 0;
}

void corpus::free_corpus()
{
    for (int i = 0; i < num_docs; i++)
    {
        document * doc = docs[i];
        delete doc;
    }
    docs.clear();
    
    size_vocab = 0;
    num_docs = 0;
    total_words = 0;
}

void corpus::read_data(const char * filename)
{
    free_corpus();
    
    int OFFSET = 0;
    FILE * fileptr;
    int length, count, word, n, nd, nw;

    // reading the data
    printf("\nreading data from %s\n", filename);

    fileptr = fopen(filename, "r");
    nd = 0;
    nw = 0;
    while ((fscanf(fileptr, "%10d", &length) != EOF))
    {
		if (length > 0) {
			document * doc = new document(length);
			for (n = 0; n < length; n++)
			{
				fscanf(fileptr, "%10d:%10d", &word, &count);
				word = word - OFFSET;
				doc->words[n] = word;
				doc->counts[n] = count;
				doc->total += count;
				if (word >= nw)
				{
					nw = word + 1;
				}
			}
			total_words += doc->total;
			doc->id = nd;
			docs.push_back(doc);
			nd++;
		}
    }
    fclose(fileptr); // close the file
    num_docs = nd;
    size_vocab = nw;
    printf("number of docs  : %d\n", nd);
    printf("number of terms : %d\n", nw);
    printf("number of total words : %d\n", total_words);
}


// end of the file
