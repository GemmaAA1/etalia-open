define(['app', 'app/model/user/user'], function (App) {

    App.Collection.Relationships = App.Backbone.Collection.extend({
        url: App.config.api_root + '/user/relationships',
        model: App.Model.Relationship
    });

    return App.Collection.Users = App.Backbone.PageableCollection.extend({
        url: App.config.api_root + '/user/users',
        model: App.Model.User
    });

});
