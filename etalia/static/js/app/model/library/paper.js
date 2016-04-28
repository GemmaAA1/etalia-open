define([
    'app',
    'app/model/library/journal',
    'app/model/library/author'
], function (App) {

    App.Model.Paper = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/library/papers',

        defaults: {
            di_doi: null,
            di_pmi: null,
            di_arx: null,
            di_pii: null,
            di_oth: null,
            title: null,
            url: null
        },

        relations: [
            {
                type: App.Backbone.HasOne,
                key: 'journal',
                relatedModel: App.Model.Journal,
                includeInJSON: 'link'
            },
            {
                type: App.Backbone.HasMany,
                key: 'authors',
                relatedModel: App.Model.Author,
                includeInJSON: 'link'
            }
        ],

        schema: {
            title: {type: 'Text', validators: ['required']},
            url: {type: 'Text', validators: ['required']}
        },

        toString: function() {
            return this.get('title');
        }
    });

    /**
     * Handlebars helpers
     */
    App.Handlebars.registerHelper('paper_title_authors', function(paper) {
        if (!paper) {
            return 'Expected paper as first argument';
        }
        var authors = paper.get('authors').map(function (author) {
            return author.get('first_name') + " " + author.get('last_name');
        });
        var output = paper.get('title') + " (" + authors.splice(0, 4).join(', ') + ")";

        return new App.Handlebars.SafeString(output);
    });

    return App.Model.Paper;
});
