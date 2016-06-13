define([
    'app',
    'app/model/library/journal',
    'app/model/library/author',
    'app/model/library/state'
], function (App) {

    var path = '/library/papers';

    App.Model.Paper = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + path,

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
                });
                this.set({state: state});
            }
            return state;
        },

        isPinned: function() {
            return this.get('state').get('watch') === App.Model.PaperState.WATCH_PINNED;
        },

        isBanned: function() {
            return this.get('state').get('watch') === App.Model.PaperState.WATCH_BANNED;
        },

        isAdded: function() {
            return this.get('state').get('store') === App.Model.PaperState.STORE_ADDED;
        },

        isTrashed: function() {
            return this.get('state').get('store') === App.Model.PaperState.STORE_TRASHED;
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
        var authors = this.authors.map(function (author) {
            return author.get('first_name') + " " + author.get('last_name');
        });
        return new App.Handlebars.SafeString(authors.join(', '));
    });
    App.Handlebars.registerHelper('paper_journal', function() {
        return new App.Handlebars.SafeString(this.journal.get('title'));
    });
    App.Handlebars.registerHelper('paper_first_seen', function() {
        return new App.Handlebars.SafeString('17 May 2016');
    });
    App.Handlebars.registerHelper('paper_new_icon', function() {
        return false ? new App.Handlebars.SafeString('<span class="new">new</span>') : ''; // TODO
    });
    App.Handlebars.registerHelper('paper_altmetric_icon', function() {
        var badge = '<span class="metric">' +
            '<span class="altmetric-embed" data-badge-type="donut" data-link-target="_blank" ' +
                'data-hide-no-mentions="true" data-badge-popover="bottom"';
        if (this.id_doi) {
            badge += ' data-doi="' + this.id_doi + '"';
        } else if (this.id_arx) {
            badge += ' data-arxiv-id="' + this.id_arx + '"';
        } else if (this.id_pmi) {
            badge += ' data-pmid="' + this.id_pmi + '"';
        } else {
            return '';
        }
        badge += '></span></span>';

        return new App.Handlebars.SafeString(badge);
    });
    App.Handlebars.registerHelper('paper_pin_class', function() {
        if (this.state && this.state.get('watch') === App.Model.PaperState.WATCH_PINNED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('paper_ban_class', function() {
        if (this.state && this.state.get('watch') === App.Model.PaperState.WATCH_BANNED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('paper_add_class', function() {
        if (this.state && this.state.get('store') === App.Model.PaperState.STORE_ADDED) {
            return ' active';
        }
        return '';
    });
    App.Handlebars.registerHelper('paper_trash_class', function() {
        if (this.state && this.state.get('store') === App.Model.PaperState.STORE_TRASHED) {
            return ' active';
        }
        return '';
    });

    return App.Model.Paper;
});
