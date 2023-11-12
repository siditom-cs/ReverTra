from transformers import LogitsWarper
import torch


class RestrictToAaLogitsWarper(LogitsWarper):
    def __init__(self, masked_input_ids: torch.LongTensor, restrict_dict: dict, filter_value: float = -float("Inf"), min_tokens_to_keep: int = 1):
        self.masked_input_ids = masked_input_ids
        self.restrict_dict = restrict_dict
        self.filter_value = filter_value
        self.min_tokens_to_keep = min_tokens_to_keep
    def __call__(self, input_ids: torch.LongTensor, scores: torch.FloatTensor) -> torch.FloatTensor:
        cur_len = input_ids.shape[-1]-1
        vocab_size = scores.shape[-1]
        if self.masked_input_ids.shape[-1] <= cur_len:
            return scores
        for bid in range(input_ids.shape[0]):
            cur_mask_input = str(int(self.masked_input_ids[bid][cur_len].item()))
            if cur_mask_input in self.restrict_dict.keys():
                restricted_words = self.restrict_dict[cur_mask_input]
                banned_indices = set(range(vocab_size))-set(self.restrict_dict[cur_mask_input])
                banned = torch.tensor([i in banned_indices for i in range(vocab_size)])
                scores[bid][banned] = self.filter_value
        return scores

