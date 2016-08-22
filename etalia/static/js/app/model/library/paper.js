define([
    'app',
    'app/model/library/journal',
    'app/model/library/author',
    'app/model/library/state'
], function (App) {

    var path = '/library/my-papers';

    App.Model.Paper = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + path,

        defaults: {
            id_doi: null,
            id_pmi: null,
            id_arx: null,
            id_pii: null,
            id_oth: null,
            title: null,
            url: null,
            new: false,
            linked_threads_count: 0
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
            },
            {
                type: App.Backbone.HasOne,
                key: 'state',
                relatedModel: App.Model.PaperState,
                includeInJSON: false,
                reverseRelation: {
                    key: 'paper',
                    type: App.Backbone.HasOne,
                    includeInJSON: 'link'
                }
            }
        ],

        schema: {
            title: {type: 'Text', validators: ['required']},
            url: {type: 'Text', validators: ['required']}
        },

        toString: function() {
            return this.get('title');
        },

        getState: function() {
            var state = this.get('state');
            if (!state) {
                state = new App.Model.PaperState({
                    user: App.getCurrentUser(),
                    paper: this
                }, {silent: true});
                this.set({state: state}, {silent: true});
            }
            return state;
        },

        isPinned: function() {
            return this.getState().get('watch') === App.Model.PaperState.WATCH_PINNED;
        },

        isBanned: function() {
            return this.getState().get('watch') === App.Model.PaperState.WATCH_BANNED;
        },

        isAdded: function() {
            return this.getState().get('store') === App.Model.PaperState.STORE_ADDED;
        },

        isTrashed: function() {
            return this.getState().get('store') === App.Model.PaperState.STORE_TRASHED;
        },

        isInLibrary: function() {
            return this.isAdded() || !this.isTrashed();
        }
    });

    App.Model.Papers = App.Backbone.Collection.extend({
        url: App.config.api_root + path,
        model: App.Model.Paper,

        initialize: function(options) {
            this.setQuery(options ? options.query : {});
        },

        setQuery: function(query) {
            this.queryParams = App._.extend({
                view: 'nested'
            }, query);

            return this;
        }
    });

    App.Model.PageablePapers = App.Backbone.PageableCollection.extend({
        url: App.config.api_root + path,
        mode: 'infinite',
        model: App.Model.Paper,

        initialize: function(options) {
            this.setQuery(options ? options.query : {});
        },

        setQuery: function(query) {
            this.queryParams = App._.extend({
                view: 'nested'
            }, query);

            return this;
        }
    });

    /**
     * Handlebars helpers
     */
    App.Handlebars.registerHelper('paper_title_authors', function(paper) {
        if (!paper) {
            throw 'Expected paper as first argument';
        }
        var authors = paper.get('authors').map(function (author) {
            return author.get('first_name') + " " + author.get('last_name');
        });
        var output = paper.get('title') + " (" + authors.splice(0, 4).join(', ') + ")";
        return new App.Handlebars.SafeString(output);
    });
    App.Handlebars.registerHelper('paper_authors', function() {
        var authors = this.hasOwnProperty('authors') ? this.authors : this.get('authors');
        authors = authors.map(function (author) {
            return author.get('first_name') + " " + author.get('last_name');
        });
        return new App.Handlebars.SafeString(authors.join(', '));
    });
    App.Handlebars.registerHelper('paper_journal', function() {
        var journal = this.hasOwnProperty('journal') ? this.journal : this.get('journal');
        if (journal) {
            return new App.Handlebars.SafeString(journal.get('title'));
        }
        return new App.Handlebars.SafeString('');
    });
    App.Handlebars.registerHelper('paper_url', function(paper) {
        if (!paper) {
            throw 'Expected paper as first argument';
        }
        return new App.Handlebars.SafeString(paper.get('url'));
    });
    App.Handlebars.registerHelper('paper_new_icon', function() {
        var is_new = this.hasOwnProperty('new') ? this.new : this.get('is_new');
        return is_new ? new App.Handlebars.SafeString('<span class="new">new</span>') : '';
    });
    App.Handlebars.registerHelper('paper_altmetric_icon', function() {
        var id_doi = this.hasOwnProperty('id_doi') ? this.id_doi : this.get('id_doi'),
            id_arx = this.hasOwnProperty('id_arx') ? this.id_arx : this.get('id_arx'),
            id_pmi = this.hasOwnProperty('id_pmi') ? this.id_pmi : this.get('id_pmi');

        var badge = '<span class="metric">' +
            '<span class="altmetric-embed" data-badge-type="donut" data-link-target="_blank" ' +
                'data-hide-no-mentions="true" data-badge-popover="bottom"';
        if (id_doi) {
            badge += ' data-doi="' + id_doi + '"';
        } else if (id_arx) {
            badge += ' data-arxiv-id="' + id_arx + '"';
        } else if (id_pmi) {
            badge += ' data-pmid="' + id_pmi + '"';
        } else {
            return '';
        }
        badge += '></span></span>';

        return new App.Handlebars.SafeString(badge);
    });
    App.Handlebars.registerHelper('paper_pin_class', function() {
        var state = this.hasOwnProperty('state') ? this.state : this.getState();
        if (state && state.get('watch') === App.Model.PaperState.WATCH_PINNED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('paper_ban_class', function() {
        var state = this.hasOwnProperty('state') ? this.state : this.getState();
        if (state && state.get('watch') === App.Model.PaperState.WATCH_BANNED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('paper_add_class', function() {
        var state = this.hasOwnProperty('state') ? this.state : this.getState();
        if (state && state.get('store') === App.Model.PaperState.STORE_ADDED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('paper_trash_class', function() {
        var state = this.hasOwnProperty('state') ? this.state : this.getState();
        if (state && state.get('store') === App.Model.PaperState.STORE_TRASHED) {
            return ' active';
        }
        return '';
    });

    return App.Model.Paper;
});
