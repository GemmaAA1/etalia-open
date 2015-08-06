import yaml
from nlp.models import Model, LSH
from library.models import Paper
from core.constants import TIME_LAPSE_CHOICES

path_to_model_file = 'scripts/build_models/models.yaml'
with open(path_to_model_file) as stream:
    MODELS = yaml.load(stream)

def build(model_name):

    check_unique_name(MODELS)
    model_args = [model for model in MODELS if model['name'] == model_name][0]

    # Initiate model
    model = Model.objects.create(**model_args)
    # dump papers data
    papers = Paper.objects.all()
    model.dump(papers)
    model.build_vocab_and_train()
    # Propagate to LSH, journalvector, papervector
    model.propagate()


def check_unique_name(models):
    names = [(i, models['name']) for i, model in enumerate(models)]
    if not len(names) == len(set(names)):
        raise AssertionError('MODELS has duplicated names')