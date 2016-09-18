define(['app', 'app/model/user/user'], function (App) {

    App.Model.Affiliation = App.Backbone.RelationalModel.extend({
        urlRoot: App.config.api_root + '/user/affiliations',

        defaults: {
            department: null,
            institution: null,
            city: null,
            country: null
        }
    });

    /**
     * Handlebars helpers.
     */
    App.Handlebars.registerHelper('user_affiliation', function(user) {
        if (!user) {
            return 'Expected user as first argument';
        }
        var affiliation = App.getProperty(user, 'affiliation');
        if (affiliation) {
            return App.getProperty(affiliation, 'institution');
        }
        return new App.Handlebars.SafeString('&nbsp;');
    });

    return App.Model.Avatar;
});
