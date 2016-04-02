define([
    'app/config',
    'backbone',
    'model/journal',
    'model/author',
    'backbone-relational'
], function (Config, Backbone, JournalModel, AuthorModel) {

    var defaults = {
        di_doi: null,
        di_pmi: null,
        di_arx: null,
        di_pii: null,
        di_oth: null,
        title: null,
        url: null
    };

    var schema = {
        title: {type: 'Text', validators: ['required']},
        url: {type: 'Text', validators: ['required']}
    };

    var PaperModel = Backbone.RelationalModel.extend({
        urlRoot: Config.api_root + '/papers',
        defaults: defaults,
        schema: schema,
        relations: [
            {
                type: Backbone.HasOne,
                key: 'journal',
                relatedModel: JournalModel,
                includeInJSON: 'link'
            },
            {
                type: Backbone.HasMany,
                key: 'authors',
                relatedModel: AuthorModel,
                includeInJSON: 'link'
            }
        ]
    });

    PaperModel.prototype.toString = function() {
        return this.get('title');
    };

    return PaperModel;
});
