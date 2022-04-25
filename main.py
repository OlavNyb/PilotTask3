import csv
import os
import shutil
import os.path
from os import path
import random

from pathlib import Path as PathlibPath

from immuneML.IO.dataset_export.AIRRExporter import AIRRExporter
from immuneML.data_model.dataset.RepertoireDataset import RepertoireDataset
from immuneML.data_model.receptor.receptor_sequence.SequenceMetadata import SequenceMetadata
from immuneML.data_model.repertoire.Repertoire import Repertoire
from immuneML.data_model.receptor.receptor_sequence.ReceptorSequence import ReceptorSequence
from immuneML.util.PathBuilder import PathBuilder
import pandas as pd


def filter_out_motif(rep, motif: str):
    return [seq for seq in rep if motif not in seq]


def simulate_motif_in_repertoires(pos, neg, motif: str, num_of_repertoires: int, repertoire_size: int,
                                  percentage_pos_in_pos_reps: float, seed):
    random.seed(seed)

    repertoires = []
    # list of boolean values that keeps track of which repertoires are positive and negative
    disease_list = []

    num_of_pos_in_pos_repertoires = percentage_pos_in_pos_reps * repertoire_size

    filtered_neg = filter_out_motif(neg, motif)

    for x in range(num_of_repertoires):
        new_rep = []
        if x < num_of_repertoires / 2:
            disease_list.append(True)
            # POSITIVE REPERTOIRE
            for i in range(int(repertoire_size)):
                if i < num_of_pos_in_pos_repertoires:
                    new_rep.append(pos[random.randint(0, len(pos) - 1)])
                else:
                    new_rep.append(filtered_neg[random.randint(0, len(filtered_neg) - 1)])
        else:
            disease_list.append(False)
            # NEGATIVE REPERTOIRE
            for i in range(int(repertoire_size)):
                new_rep.append(filtered_neg[random.randint(0, len(filtered_neg) - 1)])

        repertoires.append(new_rep)

    return repertoires, disease_list


def create_receptor_sequence(seq: str, receptor_id):
    return ReceptorSequence(amino_acid_sequence=seq,
                            identifier=receptor_id,
                            metadata=SequenceMetadata(count=1, region_type="IMGT_CDR3"))


def create_repertoire(sequence_objects, export_path, repertoire_id):
    repertoire = Repertoire.build_from_sequence_objects(sequence_objects=sequence_objects, path=export_path,
                                                        filename_base=repertoire_id,
                                                        metadata={"subject_id": repertoire_id})
    return repertoire


def main(positive_data_file_name: str, negative_data_file_name: str, given_motif: str, num_of_repertoires: int,
         repertoire_size: int,
         percentage_pos_in_pos_reps: float, output: str, seed=1):
    print("STARTED...")

    neg_seqs = []
    with open(negative_data_file_name) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            neg_seqs.append(line[0])
    file.close()

    pos_seqs = []
    with open(positive_data_file_name) as file:
        tsv_file = csv.reader(file, delimiter="\t")
        for line in tsv_file:
            pos_seqs.append(line[0])
    file.close()

    export_path = PathlibPath(os.getcwd()) / (output + "/")

    # delete path directory if already exists
    if path.exists(export_path):
        shutil.rmtree(export_path)

    PathBuilder.build(export_path)

    # create temporary directory tmp
    tmp_path = export_path / "tmp"
    PathBuilder.build(tmp_path)

    reps, disease_list = simulate_motif_in_repertoires(pos_seqs, neg_seqs, given_motif, num_of_repertoires,
                                                       repertoire_size, percentage_pos_in_pos_reps, seed)

    repertoires = [
        create_repertoire([create_receptor_sequence(seq, i) for i, seq in enumerate(r, start=1)], tmp_path, num) for
        num, r in enumerate(reps, start=1)]

    metadata_path = tmp_path / "metadata.csv"

    df = pd.DataFrame({"filename": [f"repertoire{repertoire.identifier}.npy" for repertoire in repertoires],
                       "disease": disease_list})
    df.to_csv(tmp_path / "metadata.csv", index=False)

    dataset = RepertoireDataset(repertoires=repertoires, metadata_file=metadata_path)

    AIRRExporter.export(dataset, export_path)

    # delete temporary directory tmp
    shutil.rmtree(tmp_path)

    print("FINISHED")


if __name__ == '__main__':
    positive_data_file_name, negative_data_file_name = ("SimData/few_seqs.tsv", "SimData/many_seqs.tsv")
    given_motif = "SEY"
    num_of_repertoires, repertoire_size, percentage_pos_in_pos_reps = (100, 1000, 0.5)
    output = "airr_exporter_repertoire"
    seed = 1

    main(positive_data_file_name, negative_data_file_name, given_motif, num_of_repertoires, repertoire_size,
         percentage_pos_in_pos_reps, output, seed=seed)
