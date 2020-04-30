import pysam


class Info_PARIS:
    
    def __init__(self, bam_file):
        self.blocks = [] # A block consists of two intervals representing a PARIS read
        with pysam.AlignmentFile(bam_file, 'rb') as f:
            for read in f:
                block = self.__read_to_block(read)
                if block != None:
                    self.blocks.append(block)
                
                
    def __read_to_block(self, read):
        # Only one gap is allowed
        if read.cigarstring.count('N') != 1: 
            return None
        
        # Refine read (remove indels and soft clips)
        block = read.get_blocks()
        if len(block) > 2:
            i, updated_block = 0, []
            j, cigar = 0, read.cigartuples
            while j < len(cigar):
                if cigar[j][0] == 0:
                    updated_block.append(block[i])
                    i += 1
                elif cigar[j][0] == 1 or cigar[j][0] == 2:
                    updated_block[-1] = (updated_block[-1][0], block[i][1])
                    i += 1; j += 1
                j += 1
            block = updated_block

        if len(block) != 2:
            raise ValueError('Error: Read refinement failed at: ', block)
            
        # Only retain reads with two intervals longer than 15nt
        (ll, lr), (rl, rr) = block
        if lr - ll <= 15 or rr - rl <= 15:
            return None
        
        return block
    